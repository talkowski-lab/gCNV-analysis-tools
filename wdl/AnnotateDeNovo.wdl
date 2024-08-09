version 1.0

import "Structs.wdl"

workflow AnnotateDeNovo {
  input {
    File callset   # Callset produced by the gCNV pipeline
    File pedigree  # Pedigree for the entire cohort
    File intervals # gCNV intervals

    String runtime_docker 

    Array[File] dcr_files    # denoised coverage ratio files
    Array[File] dcr_indicies # index files

    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  call DeNovo {
    input:
      callset = callset,
      pedigree = pedigree,
      intervals = intervals,

      runtime_docker = runtime_docker,

      dcr_files = dcr_files,
      dcr_indicies = dcr_indicies,

      batch_ids = batch_ids,
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

    String runtime_docker

    Array[File] dcr_files
    Array[File] dcr_indicies
  
    Array[String]? batch_ids
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = 1.5 * size(dcr_files, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 4,
    cpu_cores: 4,
    disk_gb: ceil(8.0 + input_size),
    boot_disk_gb: 16,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  Int cpus = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])

  Array[String] batches = select_first([batch_ids, []])
  File dcr_paths_file = write_lines(dcr_files)
  File batch_ids_file = if length(batches) > 0 then write_lines(batches) else '/dev/null'

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
    cat '~{write_lines(dcr_indicies)}' | xargs -- touch -c -m
    if [[ '~{batch_ids_file}' != '/dev/null' ]]; then
      paste -d '\t' '~{batch_ids_file}' '~{dcr_paths_file}' > dcrs.tsv
      dcrs=dcrs.tsv
    else
      dcrs='~{dcr_paths_file}'
    fi
    Rscript /opt/gcnv/scripts/annotate_denovo_cnv.R \
      '~{callset}' \
      '~{intervals}' \
      '~{pedigree}' \
      "${dcrs}" \
      ~{cpus} \
      'denovo_annotated_calls.bed'
    gzip 'denovo_annotated_calls.bed'
  >>>

  output {
    File denovo_annotated = "denovo_annotated_calls.bed.gz"
  }
}
