# Pipeline for counting MPRA Replicates
# Requires the same barcode_orientation as MPRAmatch
# Requires the parsed file output from MPRAmatch

workflow MPRAcount {
  Array[File] replicate_fastq #Array of replicate fastq files. Each replicate should have one file.
  Array[String] replicate_id #Identifier for each replicate listed in the replicate fastq list - should be in the same order as replicate fastqs
  Array[Pair[File,String]] fastq_id = zip(replicate_fastq, replicate_id) #Pairing the fastq file to the id
  File parsed #Output of MPRAMatch pipeline
  File acc_id
  Int? barcode_orientation = 2 #2 if MPRAMatch was run with read_a as R1 and read_b as R2, otherwise it should be 1
  Int? bc_len = 20 #If this is not changed a barcode of length 20 will be defaulted.
  String? flags = "-ECSM -A 0.05" #(Error, Cigar, Start/Stop, MD, Error Cutoff)
  String? docker_tag = "latest" #String indicating the tag of the docker image to use
  String id_out #Overall project id for the final count table
  #String working_directory #directory relative to the WDL where the scripts live
  #String out_directory #directory relative to the WDL where relevant files will be moved

  Int disk_pad = 7

  Int prep_disk = ceil(2*size(replicate_fastq[1], "GB")) + disk_pad
  Int assoc_disk = ceil(size(parsed, "GB") + 2*size(replicate_fastq[1], "GB")) + disk_pad
  #Int make_disk = ceil(length(associate.outF)*size(associate.outF[1], "GB")) + disk_pad
  Int count_disk = ceil(4*length(associate.outF)*size(associate.outF[1], "GB")) + disk_pad
  Int qc_disk = disk_pad
  Int raw_disk = disk_pad

  #Int assoc_mem = ceil(2.3*size(prep_counts.out[1], "GB"))

  scatter (replicate in fastq_id) {
    call prep_counts { input:
                          #working_directory=working_directory,
                          sample_fastq=replicate.left,
                          barcode_orientation=barcode_orientation,
                          bc_len=bc_len,
                          sample_id=replicate.right,
                          docker_tag=docker_tag,
                          prep_disk=prep_disk
                        }
    call associate { input:
                        #working_directory=working_directory,
                        matched=prep_counts.out,
                        parsed=parsed,
                        barcode_orientation=barcode_orientation,
                        sample_id=replicate.right,
                        docker_tag=docker_tag,
                        assoc_disk=assoc_disk,
                        #assoc_mem=assoc_mem
                      }
                    }
  #call make_infile { input:
  #                      #working_directory=working_directory,
  #                      tag_files=associate.outF,
  #                      tag_ids=associate.outS,
  #                      id_out=id_out,
  #                      docker_tag=docker_tag,
  #                      make_disk=make_disk
  #                    }
  call make_count_table { input:
                            #working_directory=working_directory,
                            #out_directory=out_directory,
                            tag_files=associate.outF,
                            tag_ids=associate.outS,
                            #list_inFile=make_infile.out,
                            flags=flags,
                            id_out=id_out,
                            acc_id=acc_id,
                            docker_tag=docker_tag,
                            count_disk=count_disk
                          }
  call count_QC { input:
                    count_out = make_count_table.count,
                    #out_directory = out_directory,
                    #working_directory = working_directory,
                    id_out = id_out,
                    acc_id = acc_id,
                    docker_tag=docker_tag,
                    qc_disk=qc_disk
                }
  call countRaw { input:
                    count_out = make_count_table.count,
                    cond_out = count_QC.out,
                    id_out = id_out,
                    docker_tag=docker_tag,
                    #out_directory = out_directory,
                    #working_directory = working_directory
                    raw_disk=raw_disk
                }
  # call relocate { input:
  #                   matched = prep_counts.out,
  #                   tag_files = associate.outF,
  #                   count_out = make_count_table.count,
  #                   count_log = make_count_table.log,
  #                   count_stats = make_count_table.stats,
  #                   cond_out = count_QC.out,
  #                   out_directory = out_directory
  #                 }
}

task prep_counts {
  # Grab the barcodes and check that they exist in the dictionary - if they don't exist write them to a seqparate fastq
  File sample_fastq
  Int barcode_orientation
  Int bc_len
  Int prep_disk
  #String working_directory
  String sample_id
  String docker_tag

  command {
    python /scripts/make_counts.py ${sample_fastq} ${sample_id} ${barcode_orientation} ${bc_len}
    }
  output {
    File out="${sample_id}.match"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${prep_disk} SSD"
    }
  }
task associate {
  # Associate the matched barcodes with the associated oligos
  File matched
  File parsed
  Int barcode_orientation
  Int assoc_disk
  #Int assoc_mem
  #String working_directory
  String sample_id
  String docker_tag
  command {
    perl /scripts/associate_tags.pl ${matched} ${parsed} ${sample_id}.tag ${barcode_orientation}
    }
  output {
    File outF="${sample_id}.tag"
    String outS="${sample_id}"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
    memory: "12G"
    cpu: 4
    disks: "local-disk ${assoc_disk} SSD"
    }
  }
#task make_infile {
#  # make a list of the association output and the tag ids to pass to the barcode compilation function
#  Array[File] tag_files
#  Array[String] tag_ids
#  Int make_disk
#  #String working_directory
#  String id_out
#  String docker_tag
#  command <<<
#    python /scripts/make_infile.py ${sep=',' tag_ids} ${sep=',' tag_files} ${id_out}
#  >>>
#  output {
#    File out="${id_out}_samples.txt"
#    }
#  runtime {
#    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
#    memory: "3000 MB"
#    disks: "local-disk ${make_disk} SSD"
#    }
#  }
task make_count_table {
  # Compile barcodes into a count table - columns from left to right: barcode oligo (Error CIGAR MD Aln_Start:Stop) [replicate names]
  Array[File] tag_files
  Array[String] tag_ids
  #File list_inFile
  File acc_id
  Int count_disk
  String? flags = ""
  String id_out
  String docker_tag

  Int num_rep = length(tag_ids)
  command <<<
    perl /scripts/compile_bc_cs_terra.pl ${flags} ${id_out}.count ${num_rep} ${sep=' ' tag_files} ${sep=' ' tag_ids} > ${id_out}.log
    awk '{if(NR%7==1){sum=0;good=0;bc=0;over=0;}
      if(NR%7==1){printf "%s\t",$3; printf "%s\t", ${id_out};}
      if(NR%7==3){sum+=3;bc+=$2;over+=$3;}
      if(NR%7==4){printf "%0.f\t", $2; printf "%0.f\t", $3;good+=$3;sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==5){sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==6){sum+=$3;bc+=$2;over+=$3;}
      if(NR%7==0){printf "%0.f\t", sum; printf "%.2f\t", good/(sum)*100; printf "%0.f\t", bc; printf "%0.f\t", over; printf "%.2f\n", good/(over)*100}
      }' ${id_out}.log > ${id_out}.stats
    Rscript /scripts/read_stats.R ${id_out}.stats ${acc_id} ${id_out}
    >>>
  output {
    File count="${id_out}.count"
    File log="${id_out}.log"
    File stats="${id_out}.stats"
    File plot="${id_out}_read_stats.pdf"
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
    memory: "40G"
    cpu: 10
    disks: "local-disk ${count_disk} SSD"
    }
  }
task count_QC {
  File acc_id
  File count_out
  Int qc_disk
  String id_out
  String docker_tag
  command {
    Rscript /scripts/count_QC.R ${acc_id} ${count_out} ${id_out}
    }
  output {
    File out="${id_out}_condition.txt"
    Array[File]+ plots=glob("*_QC.pdf")
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${qc_disk} SSD"
    }
  }
task countRaw {
  File count_out
  File cond_out
  Int raw_disk
  String id_out
  String docker_tag
  command {
    Rscript /scripts/bc_raw.R ${cond_out} ${count_out} ${id_out}
    }
  output {
    Array[File]+ out=glob("*.counts")
    }
  runtime {
    docker: "quay.io/tewhey-lab/mpracount:${docker_tag}"
    memory: "3000 MB"
    disks: "local-disk ${raw_disk} SSD"
    }
  }
# task relocate {
#   Array[File] matched
#   Array[File] tag_files
#   File count_out
#   File count_log
#   File count_stats
#   File cond_out
#   String out_directory
#   command <<<
#     mv ${sep=' ' matched} ${sep=' ' tag_files} ${count_out} ${count_log} ${count_stats} ${cond_out} ${out_directory}
#     >>>
#   }
