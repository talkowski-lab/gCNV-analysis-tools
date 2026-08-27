version 1.0

workflow PlotCNVEvidence {
  input {
    File callset            # GATK-gCNV callset
    File intervals          # GATK-gCNV intervals
    Array[String] batch_ids # batch IDs
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files

    String runtime_docker

    Int calls_per_split = 300 # number of calls per shard
  }

  call SplitCallset {
    input:
      callset = callset,
      calls_per_split = calls_per_split,
      runtime_docker = runtime_docker
  }

  scatter (split in SplitCallset.splits) {
    call PlotRd {
      input:
        callset = split,
        intervals = intervals,
        batch_ids = batch_ids,
        dcr_files = dcr_files,
        dcr_indices = dcr_indices,

        runtime_docker = runtime_docker
    }
  }

  call MergePlots {
    input:
      plots = PlotRd.plots,
      runtime_docker = runtime_docker
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

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Float? mem_gib
    Int? max_retries
    Int? preemptible_tries
  }

  Float inputs_size_gb = size(callset, "GB")
  Int default_disk_gb = ceil(inputs_size_gb) * 2 + 16

  Float default_mem_gib = 3.75
  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpu: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, default_disk_gb])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([mem_gib, default_mem_gib])} GiB"
    preemptible: select_first([preemptible_tries, 3])
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

  Float inputs_size_gb = size(callset, "GB") +
    size(intervals, "GB") +
    size(dcr_files, "GB") +
    size(dcr_indices, "GB")
  Int default_disk_gb = ceil(inputs_size_gb) * 2 + 16

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 2])
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

    Rscript /opt/gcnv/scripts/plot_cnv_evidence.R \
      '~{callset}' \
      '~{intervals}' \
      dcrs.tsv \
      'rd_plots'
    tar --create --gzip --file='rd_plots.tar.gz' 'rd_plots'
  >>>

  output {
    File plots = "rd_plots.tar.gz"
  }
}

task MergePlots {
  input {
    Array[File] plots
    String runtime_docker

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Float? mem_gib
    Int? max_retries
    Int? preemptible_tries
  }

  Float inputs_size_gb = size(plots, "GB")
  Int default_disk_gb = ceil(inputs_size_gb * 5) + 16
  Float default_mem_gib = 3.75

  runtime {
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    cpus: select_first([cpus, 1])
    disks: "local-disk ${select_first([disk_gb, default_disk_gb])} HDD"
    docker: runtime_docker
    maxRetries: select_first([max_retries, 1])
    memory: "${select_first([mem_gib, default_mem_gib])} GiB"
    preemptible: select_first([preemptible_tries, 3])
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
