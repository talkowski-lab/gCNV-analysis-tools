version 1.0

workflow ConcatSampleBatchMaps {
  input {
    Array[File] sample_batch_maps
    String? output_name

    String runtime_docker
  }

  call DoConcat {
    input:
      sample_batch_maps = sample_batch_maps,
      output_name = output_name,
      runtime_docker = runtime_docker
  }

  output {
    File merged_sample_batch_map = DoConcat.combined_map
  }
}

task DoConcat {
  input {
    Array[File] sample_batch_maps
    String? output_name
    String runtime_docker
  }

  Int default_disk_gb = ceil(size(sample_batch_maps, "GB") * 2) + 16

  runtime {
    bootDiskSizeGb: 8
    cpus: 1
    disks: "local-disk ~{default_disk_gb} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "3.75 GiB"
    preemptible: 3
  }

  String output_name = "~{default='combined_batch_map' output_name}.tsv"

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    cat '~{write_lines(sample_batch_maps)}' \
      | xargs cat > '~{default="combined_batch_map" output_name}.tsv'
  >>>

  output {
    File combined_map = output_name
  }
}
