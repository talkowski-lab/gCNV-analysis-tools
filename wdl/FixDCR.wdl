version 1.0

workflow FixDCR {
  input {
    File dcr
    String batch_id
    String docker
  }

  call MakeStandardDCR {
    input:
      dcr = dcr,
      batch_id = batch_id,
      runtime_docker = docker
  }

  output {
    File fixed_dcr = MakeStandardDCR.std_dcr
    File fixed_dcr_index = MakeStandardDCR.std_dcr_index
  }
}

task MakeStandardDCR {
  input {
    File dcr
    String batch_id
    String runtime_docker
  }

  Int default_disk_gb = ceil(size(dcr, "GB") * 2) + 8

  command <<<
    set -o errexit
    set -o pipefail
    set -o nounset

    mv '~{dcr}' old.dcr.bed.gz
    bgzip -dc old.dcr.bed.gz \
      | sed -e '1s/CHR/chr/' -e '1s/START/start/' -e '1s/END/end/' \
      | bgzip > '~{batch_id}.dcr.bed.gz'
    tabix -s 1 -b 2 -e 3 -c '#' '~{batch_id}.dcr.bed.gz'
  >>>

  runtime {
    bootDiskSizeGb: 8
    cpu: 1
    disks: "local-disk ${default_disk_gb} HDD"
    docker: runtime_docker
    maxRetries: 1
    memory: "3.75 GiB"
    preemptible: 3
  }

  output {
    File std_dcr = "~{batch_id}.dcr.bed.gz"
    File std_dcr_index = "~{batch_id}.dcr.bed.gz.tbi"
  }
}
