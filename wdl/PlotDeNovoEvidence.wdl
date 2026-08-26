version 1.0

import "Structs.wdl"

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

    RuntimeAttr? runtime_attr_override
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
      runtime_docker = runtime_docker,
      runtime_attr_override = runtime_attr_override
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

    RuntimeAttr? runtime_attr_override
  }

  Float input_size = (1.2 * size(dcr_files, 'GB')) + size(denovos, 'GB') + size(pedigree, 'GB') + (size(sample_batch_map, 'GB') * 2)
  Int disk_size_gb = ceil(input_size) + 16

  RuntimeAttr runtime_default = object {
    mem_gb: 8,
    cpu_cores: 2,
    disk_gb: disk_size_gb,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Int cpus = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB'
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
    File plots = 'rd_plots.tar.gz'
  }
}
