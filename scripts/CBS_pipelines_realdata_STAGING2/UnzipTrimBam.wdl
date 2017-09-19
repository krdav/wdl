
task UnzipAndSplit {
  File input_fastqR1
  File input_fastqR2
  File PIGZ
  String sample_name
  Int cpu=2

  # Uncompress while splitting fastq into chuncks of 10E7 reads:
  command {
    ${PIGZ} -dc -p 2 ${input_fastqR1} | split -l 40000000 --additional-suffix=".fastq" - "${sample_name}_1_" &
    ${PIGZ} -dc -p 2 ${input_fastqR2} | split -l 40000000 --additional-suffix=".fastq" - "${sample_name}_2_" &
    wait
  }
  runtime {
    cpu: cpu
  }
  output {
    Array[File] R1_splits = glob("*_1_??.fastq")
    Array[File] R2_splits = glob("*_2_??.fastq")
  }
}

task TrimReads {
  File input_fastqR1
  File input_fastqR2
  String basenameR1
  String basenameR2
  File TRIMMOMATIC
  File adapters = '/services/tools/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa'
  Int cpu=4

  command {
    java -Xmx80g \
      -jar ${TRIMMOMATIC} \
      PE \
      -threads ${cpu} \
      -phred33 \
      ${input_fastqR1} ${input_fastqR2} ${basenameR1}_trim.fastq ${basenameR1}_trim_unpaired.fastq ${basenameR2}_trim.fastq ${basenameR2}_trim_unpaired.fastq \
      ILLUMINACLIP:${adapters}:2:30:10
  }
  runtime {
    cpu: cpu
  }
  output {
    File output_R1 = "${basenameR1}_trim.fastq"
    File output_R2 = "${basenameR2}_trim.fastq"
  }
}

task FastqToBam {
  File input_fastqR1
  File input_fastqR2
  String basenameR1 = basename(input_fastqR1, ".fastq")
  String basename = sub(basenameR1, '_1', '')
  File PICARD
  String sample_name
  Int cpu=2

  command {
    java -Xmx8G \
      -jar ${PICARD} FastqToSam \
      FASTQ=${input_fastqR1} \
      FASTQ2=${input_fastqR2} \
      OUTPUT=${basename}.bam \
      READ_GROUP_NAME=H0164.2 \
      SAMPLE_NAME=${sample_name} \
      LIBRARY_NAME=Solexa-272222 \
      PLATFORM_UNIT=UNKNOWN \
      PLATFORM=illumina \
      SEQUENCING_CENTER=BGI \
      RUN_DATE=`date +"%y-%m-%dT%H:%M:%S"`
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${basename}.bam"
  }
}


workflow UnzipTrimBam {
  File rawdata_fastqR1
  File rawdata_fastqR2
  String sample_name
  File pigz
  File trimmomatic
  File picard

  call UnzipAndSplit {
      input:
        PIGZ=pigz,
        input_fastqR1 = rawdata_fastqR1,
        input_fastqR2 = rawdata_fastqR2,
        sample_name = sample_name
  }

  Array[Pair[File, File]] fastqR1R2_chunks = zip(UnzipAndSplit.R1_splits, UnzipAndSplit.R2_splits)
  scatter (fastq_chunk in fastqR1R2_chunks) {
    call TrimReads {
        input:
          TRIMMOMATIC=trimmomatic,
          input_fastqR1 = fastq_chunk.left,
          input_fastqR2 = fastq_chunk.right,
          basenameR1 = basename(fastq_chunk.left, ".fastq"),
          basenameR2 = basename(fastq_chunk.right, ".fastq")
    }

    call FastqToBam {
        input:
          PICARD=picard,
          sample_name = sample_name,
          input_fastqR1 = TrimReads.output_R1,
          input_fastqR2 = TrimReads.output_R2
    }
  }

  output {
    FastqToBam.out_bam,
    TrimReads.output_R1,
    TrimReads.output_R2
  }

}


