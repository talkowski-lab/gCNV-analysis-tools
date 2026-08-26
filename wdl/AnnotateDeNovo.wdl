version 1.0

import "Structs.wdl"

workflow AnnotateDeNovo {
  input {
    File sample_batch_map   # file with batch ID to sample ID mapping
    Array[String] batch_ids # batch IDs
    File callset            # Callset produced by the GATK-gCNV pipeline
    File intervals          # GATK-gCNV intervals
    File pedigree           # Pedigree for the entire cohort
    Array[File] dcr_files   # denoised copy ratio files
    Array[File] dcr_indices # denoised copy ratio index files

    String runtime_docker

    Boolean? recal_freq     # recalulate variant frequency (default is true)
    String? hq_cols         # list of columns that indicate high-quality calls
    Float? max_freq         # maximum variant frequency to consider
    Boolean? skip_allosomes # skip de novo calling on allosomes

    RuntimeAttr? runtime_attr_override
  }

  call DeNovo {
    input:
      callset = callset,
      pedigree = pedigree,
      intervals = intervals,
      sample_batch_map = sample_batch_map,
      batch_ids = batch_ids,
      dcr_files = dcr_files,
      dcr_indices = dcr_indices,
      runtime_docker = runtime_docker,
      recal_freq = recal_freq,
      hq_cols = hq_cols,
      max_freq = max_freq,
      skip_allosomes = skip_allosomes,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File denovo_annotated = DeNovo.denovo_annotated
  }
}

task DeNovo {
  input {
    File callset
    File pedigree
    File intervals
    File sample_batch_map
    Array[String] batch_ids
    Array[File] dcr_files
    Array[File] dcr_indices

    String runtime_docker

    Boolean recal_freq = true
    String? hq_cols
    Float? max_freq
    Boolean skip_allosomes = false

    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(dcr_files, 'GB') + size(dcr_indices, 'GB') + size(callset, 'GB') +
    size(pedigree, 'GB') + size(intervals, 'GB')
  Float output_size = size(callset, 'GB')
  Int disk_size_gb = ceil(input_size + output_size) + 16

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

    Rscript /opt/gcnv/scripts/annotate_denovo_cnv.R \
      ~{if recal_freq then '' else '--no-recal-freq'} \
      ~{if defined(hq_cols) then '--hq-cols ~{hq_cols}' else ''} \
      ~{if defined(max_freq) then '--max-freq ~{max_freq}' else ''} \
      ~{if skip_allosomes then '--skip-allosomes' else ''} \
      --cpus ~{cpus} \
      '~{callset}' \
      '~{intervals}' \
      '~{pedigree}' \
      dcrs.tsv \
      'denovo_annotated_calls.tsv'
    gzip 'denovo_annotated_calls.tsv'
  >>>

  output {
    File denovo_annotated = 'denovo_annotated_calls.tsv.gz'
  }
}
