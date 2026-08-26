version 1.0

import "Structs.wdl"

workflow CheckPloidy {
  input {
    File sample_batch_map   # file with batch ID to sample ID mapping
    Array[String] batch_ids # batch IDs
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files
    Array[String]? contigs  # chromosomes to check
    Boolean? is_hg19        # assume hg19 chromosome names

    String runtime_docker
  }

  Array[String] hg38_contigs = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22'
  ]
  Array[String] hg19_contigs = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',
    '14', '15', '16', '17', '18', '19', '20', '21', '22'
  ]
  Boolean hg19 = select_first([is_hg19, false])

  Array[String] default_contigs = if hg19 then hg19_contigs else hg38_contigs
  Array[String] ploidy_contigs = select_first([contigs, default_contigs])

  scatter (co in ploidy_contigs) {
    call CheckContigsPloidy {
      input:
        sample_batch_map = sample_batch_map,
        batch_ids = batch_ids,
        dcr_files = dcr_files,
        dcr_indices = dcr_indices,
        contigs = [co],
        runtime_docker = runtime_docker,
    }
  }

  call MergePloidy {
    input:
      ploidy_tars = CheckContigsPloidy.ploidy,
      runtime_docker = runtime_docker
  }

  output {
    File ploidy_plots = MergePloidy.ploidy_plots
    File ploidy_matrix = MergePloidy.ploidy_matrix
    File? aneuploidies = MergePloidy.aneuploidies
  }
}

task CheckContigsPloidy {
  input {
    File sample_batch_map
    Array[String] batch_ids
    Array[File] dcr_files
    Array[File] dcr_indices
    Array[String] contigs
    String runtime_docker

    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(read_lines(sample_batch_map))
  Float input_size = size(dcr_files, 'GB') +
    size(dcr_indices, 'GB') +
    size(sample_batch_map, 'GB')
  Float output_size = num_samples * length(contigs) * 0.000000016
  Int disk_size_gb = ceil(input_size) + ceil(output_size) + 8

  Float mem_gb = (ceil(num_samples / 1000) * 1.0) + 1.0

  RuntimeAttr runtime_default = object {
    mem_gb: mem_gb,
    cpu_cores: 1,
    disk_gb: disk_size_gb,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB'
    disks: 'local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD'
    cpus: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    dcr_paths='~{write_lines(dcr_files)}'
    cat '~{write_lines(dcr_indices)}' | xargs -- touch -c -m
    paste -d '\t' '~{write_lines(batch_ids)}' "${dcr_paths}" > dcrs.tsv

    Rscript /opt/gcnv/scripts/check_ploidy.R \
      '~{sample_batch_map}' \
      '~{write_lines(contigs)}' \
      dcrs.tsv \
      ploidy
    tar --create --file ploidy.tar ploidy
  >>>

  output {
    File ploidy = 'ploidy.tar'
  }
}

task MergePloidy {
  input {
    Array[File] ploidy_tars
    String runtime_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(ploidy_tars, 'GB')
  Float output_size = input_size * 2
  Int disk_size_gb = ceil(input_size * 2) + ceil(output_size) + 8

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    cpu_cores: 1,
    disk_gb: disk_size_gb,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB'
    disks: 'local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD'
    cpus: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    docker: runtime_docker
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mkdir ploidy
    mkdir matrices
    printf 'sample\tchr\tmean_dCR\n' > aneuploidies.tsv
    i=0
    while read -r f; do
      tar --extract --file "${f}" --wildcards '*.png'
      { tar --extract --to-stdout --file "${f}"  'ploidy/aneuploidies/aneuploidies.tsv' 2>/dev/null || true; } \
        | awk 'NR > 1' >> aneuploidies.tsv
      tar --extract --to-stdout --file "${f}" 'ploidy/ploidy_matrix.tsv' > "matrices/${i}-ploidy.tsv"
      i=$(( i + 1 ))
    done < '~{write_lines(ploidy_tars)}'

    Rscript /opt/gcnv/scripts/merge_ploidy_matrices.R matrices ploidy_matrix.tsv
    gzip ploidy_matrix.tsv
    tar --create --gzip --file ploidy_plots.tar.gz ploidy

    if [[ $(wc -l aneuploidies.tsv | awk '{print $1}') -eq 1 ]]; then
      rm aneuploidies.tsv
    fi
  >>>

  output {
    File ploidy_plots = 'ploidy_plots.tar.gz'
    File ploidy_matrix = 'ploidy_matrix.tsv.gz'
    File? aneuploidies = 'aneuploidies.tsv'
  }
}
