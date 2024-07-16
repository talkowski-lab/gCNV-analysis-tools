version 1.0

import "Structs.wdl"

workflow PlotCNVEvidence {
  input {
    File callset   # Callset produced by the gCNV pipeline or a subset of it
    File? denovo   # de novo calls
    File pedigree  # Pedigree for the entire cohort
    File intervals # gCNV intervals

    String runtime_docker

    Array[Array[File]] dcr_file_sets  # denoised coverage ratio files
    Array[Array[File]] dcr_index_sets # index files

    RuntimeAttr? runtime_attr_override
  }

  Array[File] dcr_files = flatten(dcr_file_sets)
  Array[File] dcr_indicies = flatten(dcr_index_sets)

  call PlotRD

  output {
    File plots = PlotRD.plots
  }
}

task PlotRD {
  input {
    File callset
    File? denovo
    File pedigree
    File intervals

    String runtime_docker

    Array[File] dcr_files
    Array[File] dcr_indicies

    RuntimeAttr? runtime_attr_override
  }

  Float input_size = 1.5 * size(dcr_files, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 4,
    disk_gb: ceil(8.0 + input_size),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 16
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Int cpus = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
  
  runtime {
    memory: "${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ${select_first([runtime_attr.disk_gb, runtime_default.mem_gb])} HDD"
    cpu: cpus
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: runtime_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    Rscript /opt/gcnv/scripts/plot_cnv_evidence.R \
      '~{callset}' \
      '~{denovo}' \
      '~{write_lines(dcr_files)}' \
      '~{pedigree}' \
      '~{intervals}' \
      'rd_plots'
    tar --create --gzip --file='rd_plots.tar.gz' 'rd_plots'
  >>>

  output {
    File plots = "rd_plots.tar.gz"
  }
}
