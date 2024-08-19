version 1.0

import "Structs.wdl"

workflow PlotCNVEvidence {
  input {
    File callset   # Callset produced by the gCNV pipeline or a subset of it
    File? denovo   # de novo calls
    File pedigree  # Pedigree for the entire cohort
    File intervals # gCNV intervals

    String runtime_docker

    Array[File] dcr_files    # denoised coverage ratio files
    Array[File] dcr_indicies # index files

    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  call PlotRD {
    input:
      callset = callset,
      denovo = denovo,
      pedigree = pedigree,
      intervals = intervals,

      runtime_docker = runtime_docker,

      dcr_files = dcr_files,
      dcr_indicies = dcr_indicies,

      batch_ids = batch_ids,
      runtime_attr_override = runtime_attr_override
  }

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

    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = 1.5 * size(dcr_files, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: ceil(8.0 + input_size),
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Int cpus = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])

  Array[String] batch_ids_arr = select_first([batch_ids, []])
  Boolean make_dcr_map = length(batch_ids_arr) > 0

  runtime {
    memory: "${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
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
    cat '~{write_lines(dcr_indicies)}' | xargs -- touch -c -m
    if [[ '~{make_dcr_map}' = 'true' ]]; then
      paste -d '\t' '~{write_lines(batch_ids_arr)}' "${dcr_paths}" > dcrs.tsv
      dcrs=dcrs.tsv
    else
      dcrs="${dcr_paths}"
    fi
    Rscript /opt/gcnv/scripts/plot_cnv_evidence.R \
      '~{callset}' \
      '~{denovo}' \
      '~{intervals}' \
      '~{pedigree}' \
      "${dcrs}" \
      'rd_plots'
    tar --create --gzip --file='rd_plots.tar.gz' 'rd_plots'
  >>>

  output {
    File plots = "rd_plots.tar.gz"
  }
}
