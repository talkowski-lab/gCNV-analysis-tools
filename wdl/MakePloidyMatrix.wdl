version 1.0

import "Structs.wdl"

workflow MakePloidyMatrix {
  input {
    File callset

    Array[File] dcr_files
    Array[File] dcr_indicies

    String runtime_docker

    Array[String]? contigs
    Array[String]? batch_ids

    RuntimeAttr? runtime_attr_override
  }

  call PloidyMatrix {
    input:
      callset = callset,

      dcr_files = dcr_files,
      dcr_indicies = dcr_indicies,

      runtime_docker = runtime_docker,

      contigs = contigs,
      batch_ids = batch_ids,

      runtime_attr_override = runtime_attr_override
  }

  output {
    File ploidy_matrix = PloidyMatrix.ploidy_matrix
  }
}

task PloidyMatrix {
  input {
    File callset

    Array[File] dcr_files
    Array[File] dcr_indicies

    String runtime_docker

    Array[String]? contigs
    Array[String]? batch_ids

    RuntimeAttr? runtime_attr_override
  }

  Array[String] default_contigs = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
  ]
  Array[String] ploidy_contigs = select_first([contigs, default_contigs])

  Float input_size = size(dcr_files, "GB") + size(dcr_indicies, "GB") + size(callset, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2,
    cpu_cores: 1,
    disk_gb: ceil(8.0 + input_size),
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Array[String] batch_ids_arr = select_first([batch_ids, []])
  Boolean make_dcr_map = length(batch_ids_arr) > 0

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: runtime_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    dcr_paths='~{write_lines(dcr_files)}'
    cat '~{write_lines(dcr_indicies)}' | xargs -- touch -c -m
    if [[ '~{make_dcr_map}' = 'true' ]]; then
      paste -d '\t' '~{write_lines(batch_ids_arr)}' "${dcr_paths}" > dcrs.tsv
      dcrs=dcrs.tsv
    else
      dcrs="${dcr_paths}"
    fi
    Rscript /opt/gcnv/scripts/make_ploidy_matrix.R \
      '~{callset}' \
      "${dcrs}" \
      '~{write_lines(ploidy_contigs)}' \
      'ploidy_matrix.tsv
    gzip 'ploidy_matrix.tsv'
  >>>

  output {
    File ploidy_matrix = "ploidy_matrix.tsv.gz"
  }
}
