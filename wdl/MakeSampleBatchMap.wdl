version 1.0

workflow MakeSampleBatchMap {
  input {
    Array[String] sample_ids
    Array[String] sample_set_ids
    String runtime_docker

    String? sample_set_id_trim_regex
  }

  call TrimSampleSetIDs {
    input:
      sample_set_ids = sample_set_ids,
      runtime_docker = runtime_docker,
      trim_regex = sample_set_id_trim_regex
  }

  call MakeMap {
    input:
      sample_ids = sample_ids,
      batch_ids = TrimSampleSetIDs.trimmed_ids,
      runtime_docker = runtime_docker,
  }

  output {
    Array[String] trimmed_sample_set_ids = TrimSampleSetIDs.trimmed_ids
    File sample_batch_map = MakeMap.map
  }
}

task TrimSampleSetIDs {
  input {
    Array[String] sample_set_ids
    String runtime_docker
    String trim_regex = '_((COHORT)|(CASE))$'
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
    awk -F'\t' '{sub(regex, "", $1); print $1}' \
      regex='~{trim_regex}' < '~{write_lines(sample_set_ids)}' \
      > trimmed_sample_set_ids.list
  >>>

  output {
    Array[String] trimmed_ids = read_lines('trimmed_sample_set_ids.list')
  }
}

task MakeMap {
  input {
    Array[String] sample_ids
    Array[String] batch_ids
    String runtime_docker
  }

  Array[Pair[String, String]] batch_map = zip(batch_ids, sample_ids)

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
    cat '~{write_json(batch_map)}' \
     | jq '.[] | .left as $a | .right | map([$a, .]) | .[] | @tsv' \
     > sample2batch_map.tsv
  >>>

  output {
    File map = 'sample2batch_map.tsv'
  }
}
