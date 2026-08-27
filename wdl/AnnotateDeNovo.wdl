version 1.0

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

    Boolean? recal_freq     # recalculate variant frequency (default is true)
    String? hq_cols         # list of columns that indicate high-quality calls
    Float? max_freq         # maximum variant frequency to consider
    Boolean? skip_allosomes # skip de novo calling on allosomes
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

    Int? boot_disk_gb
    Int? cpus
    Int? disk_gb
    Float? mem_gib
    Int? max_retries
    Int? preemptible_tries
  }

  Float inputs_size_gb = size(callset, "GB") +
    size(pedigree, "GB") +
    size(intervals, "GB") +
    size(sample_batch_map, "GB") +
    size(dcr_files, "GB") +
    size(dcr_indices, "GB")
  Int default_disk_gb = ceil(inputs_size_gb) + 50

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

    Rscript /opt/gcnv/scripts/annotate_denovo_cnv.R \
      ~{if recal_freq then '' else '--no-recal-freq'} \
      ~{if defined(hq_cols) then '--hq-cols ~{hq_cols}' else ''} \
      ~{if defined(max_freq) then '--max-freq ~{max_freq}' else ''} \
      ~{if skip_allosomes then '--skip-allosomes' else ''} \
      --cpus ~{runtime_cpus} \
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
