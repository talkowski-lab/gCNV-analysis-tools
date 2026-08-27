version 1.0

workflow PlotDeNovoEvidence {
  input {
    File sample_batch_map   # file with batch ID to sample ID mapping
    Array[String] batch_ids # batch IDs
    File denovos            # callset with de novos annotated
    File intervals          # GATK-gCNV intervals
    File pedigree           # pedigree for the entire cohort
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files

    String runtime_docker
  }

  call PlotRD {
    input:
      sample_batch_map = sample_batch_map,
      batch_ids = batch_ids,
      denovos = denovos,
      intervals = intervals,
      pedigree = pedigree,
      dcr_files = dcr_files,
      dcr_indices = dcr_indices,
      runtime_docker = runtime_docker
  }

  output {
    File plots = PlotRD.plots
  }
}

task PlotRD {
  input {
    File sample_batch_map
    Array[String] batch_ids
    File denovos
    File intervals
    File pedigree
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

  Float inputs_size_gb = size(sample_batch_map, "GB") +
    size(denovos, "GB") +
    size(intervals, "GB") +
    size(pedigree, "GB") +
    size(dcr_files, "GB") +
    size(dcr_indices, "GB")
  Int default_disk_gb = ceil(inputs_size_gb) * 2 + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpu: select_first([cpus, 2])
    disks: "local-disk ${select_first([disk_gb, default_disk_gb])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([mem_gib, 8])} GiB"
    preemptible: select_first([preemptible_tries, 3])
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    dcr_paths='~{write_lines(dcr_files)}'
    cat '~{write_lines(dcr_indices)}' | xargs -- touch -c -m
    paste -d '\t' '~{write_lines(batch_ids)}' "${dcr_paths}" > dcrs.tsv

    Rscript /opt/gcnv/scripts/plot_denovo_evidence.R \
      '~{denovos}' \
      '~{sample_batch_map}' \
      '~{intervals}' \
      '~{pedigree}' \
      dcrs.tsv \
      'rd_plots'
    tar --create --gzip --file='rd_plots.tar.gz' 'rd_plots'
  >>>

  output {
    File plots = "rd_plots.tar.gz"
  }
}
