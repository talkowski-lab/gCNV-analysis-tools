version 1.0

workflow CheckSex {
  input {
    Array[String] batch_ids # batch IDs
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files

    String runtime_docker
  }

  call EstimatePloidy {
    input:
      batch_ids = batch_ids,
      dcr_files = dcr_files,
      dcr_indices = dcr_indices,
      runtime_docker = runtime_docker,
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

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Float? mem_gib
    Int? max_retries
    Int? preemptible_tries
  }

  Float inputs_size_gb = size(dcr_files, "GB") + size(dcr_indices, "GB")
  Int default_disk_gb = ceil(inputs_size_gb) + 16

  Int runtime_cpus = select_first([cpus, 4])

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpu: runtime_cpus
    disks: "local-disk ${select_first([disk_gb, default_disk_gb])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([mem_gib, 16])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    dcr_paths='~{write_lines(dcr_files)}'
    cat '~{write_lines(dcr_indices)}' | xargs -- touch -c -m
    paste -d '\t' '~{write_lines(batch_ids)}' "${dcr_paths}" > dcrs.tsv

    Rscript /opt/gcnv/scripts/check_sex.R \
      --cpus ~{runtime_cpus} \
      dcrs.tsv \
      'sex_ploidy.tsv'
    gzip 'sex_ploidy.tsv'
  >>>

  output {
    File sex_ploidy = 'sex_ploidy.tsv.gz'
  }
}
