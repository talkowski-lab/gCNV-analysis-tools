version 1.0

import "Structs.wdl"

workflow PlotPloidy {
  input {
    File intervals # gCNV intervals

    String runtime_docker

    Array[File] dcr_files    # denoised coverage ratio files
    Array[File] dcr_indicies # index files

    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  call SplitIntervalsByContig {
    input:
      intervals = intervals,
      runtime_docker = runtime_docker
  }

  scatter (split in SplitIntervalsByContig.splits) {
    call PlotPloidyPerContig {
      input:
        intervals = split,
        dcr_files = dcr_files,
        dcr_indicies  = dcr_indicies,
        batch_ids = batch_ids,
        runtime_docker = runtime_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  call CombinePlots {
    input:
      plots = PlotPloidyPerContig.plots,
      runtime_docker = runtime_docker
  }

  output {
    File ploidy_plots = CombinePlots.plots
  }
}

task SplitIntervalsByContig {
  input {
    File intervals
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(intervals, 'GB') * 2) + 8

  runtime {
    memory: '128MB'
    disks: 'local-disk ${disk_size_gb} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir splits
    awk -F'\t' '{print $0 > ("splits/" $1 ".tsv")}' '~{intervals}'
  >>>

  output {
    Array[File] splits = glob('splits/*.tsv')
  }
}

task PlotPloidyPerContig {
  input {
    File intervals
    Array[File] dcr_files
    Array[File] dcr_indicies

    String runtime_docker

    Array[String]? batch_ids

    RuntimeAttr? runtime_attr_override
  }

  Int disk_size_gb = ceil(size(dcr_files, 'GiB') + size(dcr_indicies, 'GiB') + size(intervals, 'GiB')) + 8
  RuntimeAttr runtime_default = object {
    mem_gb: 2,
    cpu_cores: 1,
    disk_gb: disk_size_gb,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Array[String] batch_ids_arr = select_first([batch_ids, []])
  Boolean make_dcr_map = length(batch_ids_arr) > 0

  runtime {
    memory: '${select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GB'
    disks: 'local-disk ${select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD'
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: runtime_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    dcr_paths='~{write_lines(dcr_files)}'
    cat '~{write_lines(dcr_indicies)}' | xargs -- touch -c -m
    if [[ '~{make_dcr_map}' = 'true' ]]; then
      paste -d '\t' '~{write_lines(batch_ids_arr)}' "${dcr_paths}" > dcrs.tsv
      dcrs=dcrs.tsv
    else
      dcrs="${dcr_paths}"
    fi
    Rscript /opt/gcnv/scripts/plot_ploidy.R "${dcrs}" '~{intervals}' plots
  >>>

  output {
    Array[File] plots = glob('plots/*.png')
  }
}

task CombinePlots {
  input {
    Array[Array[File]] plots
    String runtime_docker
  }

  Array[File] flat_plots = flatten(plots)
  Int disk_size_gb = ceil(size(flat_plots, 'GB') * 2) + 8

  runtime {
    memory: '128MB'
    disks: 'local-disk ${disk_size_gb} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    mkdir ploidy_plots
    while read -r f; do
      mv "${f}" ploidy_plots
    done < '~{write_lines(flat_plots)}'
    tar -czf ploidy_plots.tar.gz ploidy_plots
  >>>

  output {
    File plots = 'ploidy_plots.tar.gz'
  }
}
