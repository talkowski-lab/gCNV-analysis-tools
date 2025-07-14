version 1.0

workflow MakeSampleBatchMap {
  input {
    String sample_set_id
    Array[String] sample_ids
    String runtime_docker

    String sample_set_id_trim_regex = "_((COHORT)|(CASE))$"
  }

  call MakeMap {
    input:
      sample_set_id = sample_set_id,
      sample_ids = sample_ids,
      runtime_docker = runtime_docker
  }

  output {
    File map = MakeMap.map
  }
}

task MakeMap {
  input {
    String sample_set_id
    Array[String] sample_ids
    String sample_set_id_trim_regex
    String runtime_docker
  }

  runtime {
    memory: '1 GB'
    disks: 'local-disk 16 HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    batch_id='~{sample_set_id}'
    if [[ -n '~{sample_set_id_trim_regex}' ]]; then
      batch_id="$(awk -vbatch_id="${batch_id}" -vregex="${sample_set_id_trim_regex}" 'BEGIN{sub(regex, "", batch_id); print batch_id}')"
    fi
    awk '{print batch_id "\t" $1}' batch_id="${batch_id}" '~{write_lines(sample_ids)}' > sample_batch_map.tsv
  >>>

  output {
    File map = 'sample_batch_map.tsv'
  }
}
