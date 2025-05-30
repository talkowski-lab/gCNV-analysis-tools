version 1.0

import "Structs.wdl"

workflow MergeDcr {
  input {
    Array[File] sample_dcrs
    String sample_set_id
    String runtime_docker
    RuntimeAttr? runtime_override_merge
  }

  call Merge {
    input:
      sample_dcrs = sample_dcrs,
      sample_set_id = sample_set_id,
      runtime_docker = runtime_docker,
      runtime_attr_override = runtime_override_merge
  }

  output {
    File dcr = Merge.dcr
    File dcr_index = Merge.dcr_index
  }
}

task Merge {
  input {
    Array[File] sample_dcrs
    String sample_set_id
    String runtime_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(sample_dcrs, 'GB')
  Int extra_mem = if length(sample_dcrs) <= 300 then 0 else floor(length(sample_dcrs) / 300) - 1
  RuntimeAttr runtime_default = object {
    mem_gb: 16 + extra_mem,
    cpu_cores: 2,
    disk_gb: ceil(input_size * 8) + 16,
    boot_disk_gb: 8,
    preemptible_tries: 2,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB'
    disks: 'local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD'
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: runtime_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  String output_name = sample_set_id + "_COHORT.dcr.bed.gz"
  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir dcrs
    while read -r f; do
      mv "${f}" dcrs
    done < '~{write_lines(sample_dcrs)}'

cat > commands.sql <<EOF
CREATE TABLE dcr_mat AS
SELECT * FROM read_csv('dcrs/*.tsv',
                       comment = '@',
                       delim = '\t',
                       skip = 1,
                       filename = true,
                       columns = {'chr': 'VARCHAR',
                                  'start': 'UINTEGER',
                                  'end': 'UINTEGER',
                                  'lcr': 'DOUBLE'});
COPY (
    SELECT filename, count(*) AS records FROM dcr_mat GROUP BY filename
) TO 'counts.tsv' (FORMAT csv, DELIMITER '\t', HEADER false);
EOF

    duckdb -bail dcr.duckdb < commands.sql
    declare -i prev count fnr=1
    while IFS=$'\t' read -r filename count; do
      if (( fnr == 1 )); then
        prev="${count}"
        fnr+=1
        continue
      fi

      if (( count != prev )); then
        printf '%s does not have same number of records as others\n' "${filename}" >&2
        exit 1
      fi
      fnr+=1
    done < counts.tsv

cat > commands.sql <<EOF
ALTER TABLE dcr_mat ADD COLUMN sample_id VARCHAR;
UPDATE dcr_mat
SET sample_id = regexp_replace(parse_filename(filename, true), '^denoised_copy_ratios-', '');
ALTER TABLE dcr_mat DROP filename;
CREATE TABLE dcr_mat_pivot AS
SELECT * FROM (
    PIVOT dcr_mat ON sample_id USING first(lcr) ORDER BY chr, start
) pivot_alias;
COPY (SELECT count(*) FROM dcr_mat_pivot)
TO 'all_count' (FORMAT csv, HEADER false);
COPY (SELECT count(*) FROM dcr_mat_pivot WHERE COLUMNS(*) IS NOT NULL)
TO 'not_null_count' (FORMAT csv, HEADER false);
COPY dcr_mat_pivot TO '/dev/stdout' (FORMAT csv, DELIMITER '\t', HEADER true);
EOF

    duckdb -bail dcr.duckdb < commands.sql \
      | awk 'FNR==1{$0 = "#" $0} 1' \
      | bgzip -c > '~{output_name}'
    tabix --sequence 1 --begin 2 --end 3 --comment '#' --zero-based '~{output_name}'

    declare -i all_count not_null_count
    read -r all_count < all_count
    read -r not_null_count < not_null_count
    if (( all_count != not_null_count )); then
      printf 'dCR matrices have mismatched intervals\n' >&2
      exit 1
    fi

    rm dcr.duckdb
  >>>

  output {
    File dcr = output_name
    File dcr_index = output_name + ".tbi"
    File? database = "dcr.duckdb"
  }
}
