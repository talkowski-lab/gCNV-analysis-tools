version 1.0

import "Structs.wdl"

workflow CheckSex {
  input {
    Array[String] batch_ids # batch IDs
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files

    String runtime_docker

    RuntimeAttr? runtime_attr_override
  }

  call EstimatePloidy {
    input:
      batch_ids = batch_ids,
      dcr_files = dcr_files,
      dcr_indices = dcr_indices,
      runtime_docker = runtime_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File sex_ploidy = EstimatePloidy.sex_ploidy
  }
}

task EstimatePloidy {
  input {
    Array[String] batch_ids
    Array[File] dcr_files
    Array[File] dcr_indices

    String runtime_docker

    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dcr_files, 'GB') + size(dcr_indices, 'GB')
  Int disk_size_gb = ceil(input_size) + 16

  RuntimeAttr runtime_default = object {
    mem_gb: 16,
    cpu_cores: 4,
    disk_gb: disk_size_gb,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Int cpus = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB'
    disks: 'local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD'
    cpu: cpus
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
    cat '~{write_lines(dcr_indices)}' | xargs -- touch -c -m
    paste -d '\t' '~{write_lines(batch_ids)}' "${dcr_paths}" > dcrs.tsv

    Rscript /opt/gcnv/scripts/check_sex.R \
      --cpus ~{cpus} \
      dcrs.tsv \
      'sex_ploidy.tsv'
    gzip 'sex_ploidy.tsv'
  >>>

  output {
    File sex_ploidy = 'sex_ploidy.tsv.gz'
  }
}
