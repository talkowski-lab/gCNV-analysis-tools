version 1.0

workflow MakeSampleBatchMap {
  input {
    File callset
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
      callset = callset,
      runtime_docker = runtime_docker
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
    File callset
    String runtime_docker
  }

  Int disk_size_gb = ceil(size(callset, 'GB')) + 8

  runtime {
    memory: '1 GB'
    disks: 'local-disk ${disk_size_gb} HDD'
    cpus: 1
    preemptible: 3
    maxRetries: 1
    docker: runtime_docker
    bootDiskSizeGb: 16
  }

  command <<<
    make_table() {
      awk -F'\t' 'NR==1{for(i=1;i<=NF;++i){a[$i]=i}} NR>1{print $(a["batch"])"\t"$(a["sample"])}' - | LC_ALL=C sort -u
    }

    if [[ '~{callset}' = *.gz ]]; then
      gzip -cd '~{callset}' | make_table > sample2batch_map.tsv
    else
      cat '~{callset}' | make_table > sample2batch_map.tsv
    fi
  >>>

  output {
    File map = 'sample2batch_map.tsv'
  }
}
