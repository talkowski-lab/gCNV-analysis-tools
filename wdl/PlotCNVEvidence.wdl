version 1.0

import "Structs.wdl"

workflow PlotCNVEvidence {
  input {
    File callset
    Int calls_per_split = 300

    File intervals # gCNV intervals

    String runtime_docker

    Array[File] dcr_files    # denoised coverage ratio files
    Array[File] dcr_indicies # index files

    Array[String]? batch_ids
    RuntimeAttr? runtime_override_split_callset
    RuntimeAttr? runtime_override_plot_rd
    RuntimeAttr? runtime_override_merge_plots
  }

  call SplitCallset {
    input:
      callset = callset,
      calls_per_split = calls_per_split,
      runtime_docker = runtime_docker,
      runtime_attr_override = runtime_override_split_callset
  }

  scatter (split in SplitCallset.splits) {
    call PlotRd {
      input:
        callset = split,
        intervals = intervals,

        runtime_docker = runtime_docker,

        dcr_files = dcr_files,
        dcr_indicies = dcr_indicies,

        batch_ids = batch_ids,
        runtime_attr_override = runtime_override_plot_rd
    }
  }

  call MergePlots {
    input:
      plots = PlotRd.plots,
      runtime_docker = runtime_docker,
      runtime_attr_override = runtime_override_merge_plots
  }

  output {
    File plots = MergePlots.plots_tar
  }
}

task SplitCallset {
  input {
    File callset
    Int calls_per_split
    String runtime_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(callset, "GB")

  RuntimeAttr runtime_default = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 2) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
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

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    if (( ~{calls_per_split} <= 0 )); then
      printf 'calls per split must be a positive integer\n' >&2
      exit 1
    fi

    mkdir splits
    # assume callset has a header, which it must
    awk '
      NR == 1 {
        header = $0
        i = 0
        next
      }
      (NR - 1) % max == 1 {
        if (out) {
          close(out)
        }
        out = sprintf("splits/%06d.tsv", i++)
        print header > out
      }
      {
        print > out
      }' max=~{calls_per_split} '~{callset}'
  >>>

  output {
    Array[File] splits = glob("splits/*.tsv")
  }
}

task PlotRd {
  input {
    File callset
    File intervals

    String runtime_docker

    Array[File] dcr_files
    Array[File] dcr_indicies

    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dcr_files, 'GB') + size(callset, 'GB')
  Int disk_size_gb = ceil(input_size) * 2 + 16

  RuntimeAttr runtime_default = object {
    mem_gb: 4,
    cpu_cores: 1,
    disk_gb: disk_size_gb,
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 0
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
      '~{intervals}' \
      "${dcrs}" \
      'rd_plots'
    tar --create --gzip --file='rd_plots.tar.gz' 'rd_plots'
  >>>

  output {
    File plots = 'rd_plots.tar.gz'
  }
}

task MergePlots {
  input {
    Array[File] plots
    String runtime_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(plots, "GB")

  RuntimeAttr runtime_default = object {
    mem_gb: 1,
    cpu_cores: 1,
    disk_gb: ceil(input_size * 5) + 16,
    boot_disk_gb: 16,
    preemptible_tries: 3,
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

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail
    shopt -s failglob

    mkdir rd_plots
    while read -r f; do
      dir="$(dirname "${f}")"
      tar --extract --gzip --directory "${dir}" --file="${f}"
      cp "${dir}/rd_plots/"*.png rd_plots
      rm -r "${dir}/rd_plots"
    done < '~{write_lines(plots)}'

    tar --create --gzip --file='rd_plots.tar.gz' rd_plots
  >>>

  output {
    File plots_tar = "rd_plots.tar.gz"
  }
}
