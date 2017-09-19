## Copyright Broad Institute, 2017
## 
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF 
## generation) according to the GATK Best Practices (June 2016) for germline SNP and 
## Indel discovery in human whole-genome sequencing (WGS) data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# TASK DEFINITIONS

task STAR_Map {
  File STAR
  String STARindexDir
  String sample_name
  File input_fastqR1
  File input_fastqR2
#  String suffix="Aligned.sortedByCoord.out"
  String suffix="Aligned.out"
  Int cpu=28

  command {
    ln -s ${STARindexDir} GenomeDir
    # Use default mode --genomeDir=./GenomeDir/
    ${STAR} --readFilesIn ${input_fastqR1} ${input_fastqR2} \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFileNamePrefix ${sample_name} \
        --twopassMode Basic \
        --runThreadN ${cpu}
  }
  output {
    File out_bam = "${sample_name}${suffix}.bam"
    String out_sample_name = "${sample_name}${suffix}"
  }
  runtime {
    cpu: cpu
  }
}

task SplitNCigarReads {
  File GATK
  String sample_name
  String ref_fasta
  File in_bam
  String suffix="_splitNcigar"
  Int cpu=2

  command {
    # Run options are set to follow the GATK best practice for RNAseq. data.
    java -jar ${GATK} \
      -T SplitNCigarReads \
      -R ${ref_fasta}.fa \
      -I ${in_bam} \
      -o ${sample_name}${suffix}.bam \
      -rf ReassignOneMappingQuality \
      -RMQF 255 \
      -RMQT 60 \
      -U ALLOW_N_CIGAR_READS
 }
  output {
    File out_bam = "${sample_name}${suffix}.bam"
    File out_bai = "${sample_name}${suffix}.bai"
    String out_sample_name = "${sample_name}${suffix}"
  }
  runtime {
    cpu: cpu
  }
}


### This was written for the old RNAseq pipeline, but should be deprecated:
task BaseRecalibrator_RNA {
  File GATK
  String sample_name
  String ref_fasta
  File in_bam
  File in_bai
  File dbsnp
  String suffix="_recal"
  Int cpu=28

  command {
    java -jar ${GATK} \
      -T BaseRecalibrator \
      -R ${ref_fasta}.fa \
      -I ${in_bam} \
      -knownSites ${dbsnp} \
      -o ${sample_name}${suffix}.grp \
      -nct ${cpu}
 }
  output {
    File out_grp = "${sample_name}${suffix}.grp"
    String out_sample_name = "${sample_name}${suffix}"
  }
  runtime {
    cpu: cpu
  }
}

task VariantFiltration_RNA {
  File GATK
  String sample_name
  String ref_fasta
  File in_vcf
  String suffix="_filter"
  Int cpu=2

  command {
    java -jar ${GATK} \
      -T VariantFiltration \
      -R ${ref_fasta}.fa \
      -V ${in_vcf} \
      -window 35 \
      -cluster 3 \
      -filterName FS \
      -filter "FS > 3.0" \
      -filterName QD \
      -filter "QD < 2.0" \
      -o ${sample_name}${suffix}.vcf
 }
  output {
    File out_vcf = "${sample_name}${suffix}.vcf"
    String out_sample_name = "${sample_name}${suffix}"
  }
  runtime {
    cpu: cpu
  }
}

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

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  File in_bam
  String metrics_filename
  Int cpu=1
  File PICARD

  command {
    java -Xmx128m \
      -jar ${PICARD} \
      CollectQualityYieldMetrics \
      INPUT=${in_bam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  }
  runtime {
    cpu: cpu
  }
  output {
    File metrics = "${metrics_filename}"
  }
}


## Check the assumption that the final GVCF filename that is going to be used ends with .g.vcf.gz 
#task CheckFinalVcfExtension {
#  String vcf_filename
#  Int cpu=1
#  File PYTHON2
#
#  command <<<
#    ${PYTHON2} <<CODE
#    import os
#    import sys
#    filename="${vcf_filename}"
#    if not filename.endswith(".g.vcf.gz"):
#      raise Exception("input","output GVCF filename must end with '.g.vcf.gz', found %s"%(filename))
#      sys.exit(1) 
#    CODE
#  >>>
#  runtime {
#    cpu: cpu
#  }
#  output {
#    String common_suffix=read_string(stdout())
#  }
#}


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File in_bam
  String bwa_commandline
  String out_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit), 
  # listing the reference contigs that are "alternative". 
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  Int cpu=28
  File PICARD
  File SAMTOOLS

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it, we proceed
    if [ -s ${ref_alt} ]; then
      java -Xmx3000m \
        -jar ${PICARD} \
        SamToFastq \
        INPUT=${in_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ${bwa_commandline} /dev/stdin -  2> >(tee ${out_bam_basename}.bwa.stderr.log >&2) | \
      ${SAMTOOLS} view -1 - > ${out_bam_basename}.bam

      grep -m1 "read .* ALT contigs" ${out_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found, so we bail out
    else
      exit 1;
    fi
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_bam_basename}.bam"
    File bwa_stderr_log = "${out_bam_basename}.bwa.stderr.log"
  }
}

# Merge original input uBAM file with BWA-aligned BAM file
task MergeBamAlignment {
  File unmapped_bam
  File aligned_bam
  String out_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Int cpu=1
  File PICARD

  command {
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    java -Xmx2500m \
      -jar ${PICARD} \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${out_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_bam_basename}.bam"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  File in_bam
  String out_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD

  command {
    set -o pipefail

    java -Xmx4000m \
      -jar ${PICARD} \
      SortSam \
      INPUT=${in_bam} \
      OUTPUT=/dev/stdout \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false | \
    java -Xmx500m \
      -jar ${PICARD} \
      SetNmAndUqTags \
      INPUT=/dev/stdin \
      OUTPUT=${out_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      REFERENCE_SEQUENCE=${ref_fasta}
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_bam_basename}.bam"
    File out_bam_index = "${out_bam_basename}.bai"
    File out_bam_md5 = "${out_bam_basename}.bam.md5"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  File in_bam
  String out_bam_prefix
  Int cpu=1
  File PICARD

  command {
    java -Xmx5000m \
      -jar ${PICARD} \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      OUTPUT=${out_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectBaseDistributionByCycle" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="MeanQualityByCycle" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="ALL_READS"

    touch ${out_bam_prefix}.insert_size_metrics
    touch ${out_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    cpu: cpu
  }
  output {
    File base_distribution_by_cycle_pdf = "${out_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "${out_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "${out_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${out_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "${out_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "${out_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "${out_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${out_bam_prefix}.quality_distribution_metrics"
  }
}

# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
  File in_bam
  File in_bam_index
  String out_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD

  command {
    java -Xmx5000m \
      -jar ${PICARD} \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${out_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="READ_GROUP"
  }
  runtime {
    cpu: cpu
  }
  output {
    File alignment_summary_metrics = "${out_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "${out_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${out_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${out_bam_prefix}.gc_bias.summary_metrics"
  }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
  File in_bam
  File in_bam_index
  String out_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD

  command {
    java -Xmx5000m \
      -jar ${PICARD} \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=${out_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectAlignmentSummaryMetrics" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="CollectSequencingArtifactMetrics" \
      PROGRAM="CollectGcBiasMetrics" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="SAMPLE" \
      METRIC_ACCUMULATION_LEVEL="LIBRARY"

    touch ${out_bam_prefix}.insert_size_metrics
    touch ${out_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    cpu: cpu
  }
  output {
    File alignment_summary_metrics = "${out_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "${out_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "${out_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "${out_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "${out_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${out_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "${out_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${out_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "${out_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "${out_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "${out_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "${out_bam_prefix}.quality_distribution_metrics"
  }
}

# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  Array[File] in_bams
  Array[File] in_bam_indexes
  File haplotype_database_file # if this file is empty (0-length) the workflow should not do fingerprint comparison (as there are no fingerprints for the sample)
  String metrics_filename
  Int cpu=1
  File PICARD

  command <<<
    if [ -s ${haplotype_database_file} ]; then
      java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx2000m \
       -jar ${PICARD} \
       CrosscheckReadGroupFingerprints \
       OUTPUT=${metrics_filename} \
       HAPLOTYPE_MAP=${haplotype_database_file} \
       EXPECT_ALL_READ_GROUPS_TO_MATCH=true \
       INPUT=${sep=' INPUT=' in_bams} \
       LOD_THRESHOLD=-20.0
    else
      echo "No haplotype_database_file. Skipping Fingerprint check."
      touch ${metrics_filename}
    fi
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

# Check that the fingerprint of the sample BAM matches the sample array
task CheckFingerprint {
  File in_bam
  File in_bam_index
  File haplotype_database_file
  File genotypes
  String output_basename
  String sample
  Int cpu=1
  File PICARD

  command <<<
  if [ -s ${genotypes} ]; then
    java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx1024m  \
    -jar ${PICARD} \
    CheckFingerprint \
    INPUT=${in_bam} \
    OUTPUT=${output_basename} \
    GENOTYPES=${genotypes} \
    HAPLOTYPE_MAP=${haplotype_database_file} \
    SAMPLE_ALIAS="${sample}" \
    IGNORE_READ_GROUPS=true
  else
    echo "No fingerprint found. Skipping Fingerprint check."
    # We touch the outputs here in order to create 0 length files.  
    # Otherwise the task will fail since the expected outputs are not to be found.
    touch ${output_basename}.fingerprinting_summary_metrics
    touch ${output_basename}.fingerprinting_detail_metrics
  fi
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File summary_metrics = "${output_basename}.fingerprinting_summary_metrics"
    File detail_metrics = "${output_basename}.fingerprinting_detail_metrics"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  Array[File] in_bams
  String out_bam_basename
  String metrics_filename
  Int cpu=1
  File PICARD

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Xmx4000m \
      -jar ${PICARD} \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' in_bams} \
      OUTPUT=${out_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"
      CREATE_MD5_FILE=true
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  File ref_dict
  Int cpu=1

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    cpu: cpu
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  File in_bam
  File in_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File GATK

  command {
    rand=`shuf -i 1-10000000 -n 1`
    mv ${write_lines(sequence_group_interval)} $rand.intervals

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx4000m \
      -jar ${GATK} \
      -T BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${in_bam} \
      --useOriginalQualities \
      -o ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${sep=" -knownSites " known_indels_sites_VCFs} \
      -L $rand.intervals
  }
  runtime {
    cpu: cpu
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  File in_bam
  File in_bam_index
  String out_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File GATK4

  command {
    rand=`shuf -i 1-10000000 -n 1`
    mv ${write_lines(sequence_group_interval)} $rand.intervals

    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
      -jar ${GATK4} ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${in_bam} \
      --useOriginalQualities \
      -O ${out_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 \
      -L $rand.intervals
  }
  runtime {
    cpu: cpu
  }
  output {
    File recalibrated_bam = "${out_bam_basename}.bam"
    File recalibrated_bam_checksum = "${out_bam_basename}.bam.md5"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int cpu=1
  File GATK

  command {
    java -Xmx3000m \
      -cp ${GATK} org.broadinstitute.gatk.tools.GatherBqsrReports \
      I=${sep=' I=' input_bqsr_reports} \
      O=${output_report_filename}
    }
  runtime {
    cpu: cpu
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  Array[File] in_bams
  String out_bam_basename

  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' in_bams} \
      OUTPUT=${out_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_bam_basename}.bam"
    File out_bam_index = "${out_bam_basename}.bai"
    File out_bam_md5 = "${out_bam_basename}.bam.md5"
  }
}

# Validate the output bam file
task ValidateSamFile {
  File in_bam
  File in_bam_index
  String report_filename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int? max_output
  Array[String]? ignore
  Int cpu=1
  File PICARD

  command {
    java -Xmx4000m \
      -jar ${PICARD} \
      ValidateSamFile \
      INPUT=${in_bam} \
      OUTPUT=${report_filename} \
      REFERENCE_SEQUENCE=${ref_fasta} \
      ${"MAX_OUTPUT=" + max_output} \
      IGNORE=${default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
    cpu: cpu
  }
  output {
    File report = "${report_filename}"
  }
}

# Collect WGS metrics (Broad Genomics stringent QC thresholds)
task CollectWgsMetrics {
  File in_bam
  File in_bam_index
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectWgsMetrics \
      INPUT=${in_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      INTERVALS=${wgs_coverage_interval_list} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    cpu: cpu
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
  File in_bam
  File in_bam_index
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectRawWgsMetrics \
      INPUT=${in_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      INTERVALS=${wgs_coverage_interval_list} \
      OUTPUT=${metrics_filename}
  }
  runtime {
    cpu: cpu
  }
  output {
    File metrics = "${metrics_filename}"
  }
}

# Generate a checksum per readgroup
task CalculateReadGroupChecksum {
  File in_bam
  File in_bam_index
  String read_group_md5_filename
  Int cpu=1
  File PICARD

  command {
    java -Xmx1000m \
      -jar ${PICARD} \
      CalculateReadGroupChecksum \
      INPUT=${in_bam} \
      OUTPUT=${read_group_md5_filename}
  }
  runtime {
    cpu: cpu
  }
  output {
    File md5_file = "${read_group_md5_filename}"
  }
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  File in_bam
  File in_bam_index
  File contamination_sites_vcf
  File contamination_sites_vcf_index
  String output_prefix
  Int cpu=1
  File verifyBamID

  # Having to do this as a 2-step command in heredoc syntax, adding a python step to read the metrics
  # This is a hack until read_object() is supported by Cromwell.
  # It relies on knowing that there is only one data row in the 2-row selfSM TSV file
  # Piping output of verifyBamId to /dev/null so only stdout is from the python command
  command <<<
    set -e

    ${verifyBamID} \
    --verbose \
    --ignoreRG \
    --vcf ${contamination_sites_vcf} \
    --out ${output_prefix} \
    --bam ${in_bam} \
    1>/dev/null

    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
    # A zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event. 
    # If the bam isn't really empty, this is probably due to the use of a incompatible reference build between 
    # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/0.75)
        i = i + 1
    # There should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
    # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File selfSM = "${output_prefix}.selfSM"
    File depthSM = "${output_prefix}.depthSM"
    File log = "${output_prefix}.log"

    # We would like to do the following, however:
    # The object is read as a string
    # explicit string->float coercion via float(), as shown below, is supported by Cromwell
    # the interim value cannot be stored as a string and then assigned to a float. Variables intialized in output cannot be dereferenced in output.
    # Float contamination = float(read_object(${output_prefix} + ".selfSM").FREEMIX) / 0.75

    # In the interim, get the value from the python hack above:
    Float contamination = read_float(stdout())
  }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File in_bam
  File in_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int cpu=1
  File GATK

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar ${GATK} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${in_bam} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead \
      -L ${interval_list}
  }
  runtime {
    cpu: cpu
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}


task HaplotypeCaller_RNA {
  File in_bam
  File in_bam_index
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int cpu=1
  File GATK

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar ${GATK} \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${in_bam} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -stand_call_conf 20 \
      -dontUseSoftClippedBases \
      -contamination ${default=0 contamination} \
      --read_filter OverclippedRead \
      -L ${interval_list}
 }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
  runtime {
    cpu: cpu
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int cpu=1
  File PICARD

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xmx2g \
      -jar ${PICARD} \
      MergeVcfs \
      INPUT=${sep=' INPUT=' input_vcfs} \
      OUTPUT=${output_vcf_name}
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
  runtime {
    cpu: cpu
  }
}

# Validate a GVCF with -gvcf specific validation
task ValidateGVCF {
  File input_vcf
  File input_vcf_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_index
  File wgs_calling_interval_list
  Int cpu=1
  File GATK

  command {
    java -Xmx8g \
      -jar ${GATK} \
      -T ValidateVariants \
      -V ${input_vcf} \
      -R ${ref_fasta} \
      -gvcf \
      --validationTypeToExclude ALLELES \
      --reference_window_stop 208 -U  \
      --dbsnp ${dbSNP_vcf} \
      -L ${wgs_calling_interval_list}
  }
  runtime {
    cpu: cpu
  }
}

# Collect variant calling metrics from GVCF output
task CollectGvcfCallingMetrics {
  File input_vcf
  File input_vcf_index
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_index
  File ref_dict
  File wgs_evaluation_interval_list
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectVariantCallingMetrics \
      INPUT=${input_vcf} \
      OUTPUT=${metrics_basename} \
      DBSNP=${dbSNP_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
      TARGET_INTERVALS=${wgs_evaluation_interval_list} \
      GVCF_INPUT=true
  }
  runtime {
    cpu: cpu
  }
  output {
    File summary_metrics = "${metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "${metrics_basename}.variant_calling_detail_metrics"
  }
}

# Convert BAM file to CRAM format for validation
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
  File in_bam
  File ref_fasta
  File ref_fasta_index
  String output_basename
  Int cpu=1
  File SAMTOOLS
  File seq_cache_populate

  command <<<
    set -e
    set -o pipefail

    ${SAMTOOLS} view -C -T ${ref_fasta} ${in_bam} | \
    tee ${output_basename}.cram | \
    md5sum | awk '{print $1}' > ${output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    ${seq_cache_populate} -root ./ref/cache ${ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    ${SAMTOOLS} index ${output_basename}.cram
    mv ${output_basename}.cram.crai ${output_basename}.crai
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File output_cram = "${output_basename}.cram"
    File output_cram_index = "${output_basename}.crai"
    File output_cram_md5 = "${output_basename}.cram.md5"
  }
}

# Convert a CRAM file to BAM format
task CramToBam {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File cram_file
  String output_basename
  Int cpu=1
  File SAMTOOLS

command <<<
  set -e
  set -o pipefail

  ${SAMTOOLS} view -h -T ${ref_fasta} ${cram_file} |
  ${SAMTOOLS} view -b -o ${output_basename}.bam -
  ${SAMTOOLS} index -b ${output_basename}.bam
  mv ${output_basename}.bam.bai ${output_basename}.bai 
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${output_basename}.bam"
    File out_bam_index = "${output_basename}.bai"
  }
}

# WORKFLOW DEFINITION
workflow PairedEndSingleSampleWorkflow {

  File contamination_sites_vcf
  File contamination_sites_vcf_index
  File fingerprint_genotypes_file # if this file is empty (0-length) the workflow should not do fingerprint comparison (as there are no fingerprints for the sample)
  File haplotype_database_file  
  File wgs_evaluation_interval_list
  File wgs_coverage_interval_list
  
  String sample_name
  String base_file_name_normal
  String base_file_name_tumor
  String final_gvcf_ext
#  Array[File] flowcell_unmapped_bams
  File rawdata_normal_fastqR1
  File rawdata_normal_fastqR2
  File rawdata_tumor_fastqR1
  File rawdata_tumor_fastqR2
  String unmapped_bam_suffix

  Array[File] scattered_calling_intervals
  File wgs_calling_interval_list
  
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac
  
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  String recalibrated_bam_basename = base_file_name_normal + ".aligned.duplicates_marked.recalibrated"

  # Tools:
  File picard
  File gatk
  File gatk4
  File python2
  File python3
  File samtools
  File bwa
  File verifyBamID
  File seq_cache_populate
  File pigz
  File trimmomatic
  File star

  String bwa_commandline = bwa + " mem -K 100000000 -p -v 3 -t 28 -Y $bash_ref_fasta"

  String final_gvcf_name_normal = base_file_name_normal + final_gvcf_ext


  # Check that the GVCF output name follows convention
#  call CheckFinalVcfExtension {
#    input:
#      PYTHON2=python2,
#      vcf_filename = final_gvcf_name_normal
#   }


  ###### Here the workflow is adapted to .fastq files:
  # Decompress and split the files into chunks, return an array of .fastq files:
  call UnzipAndSplit {
      input:
        PIGZ=pigz,
        input_fastqR1 = rawdata_normal_fastqR1,
        input_fastqR2 = rawdata_normal_fastqR2,
        sample_name = sample_name + '_normal'
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
          sample_name = sample_name + '_normal',
          input_fastqR1 = TrimReads.output_R1,
          input_fastqR2 = TrimReads.output_R2
    }

  }

  ######


  # Align flowcell-level unmapped input bams in parallel
#  scatter (unmapped_bam in flowcell_unmapped_bams) {
  scatter (unmapped_bam in FastqToBam.out_bam) {  
    # Because of a wdl/cromwell bug this is not currently valid so we have to sub(sub()) in each task
    # String base_name = sub(sub(unmapped_bam, "gs://.*/", ""), unmapped_bam_suffix + "$", "")

    String sub_strip_path = "/home/projects/dp_00005/data/cromwell_test/.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"

    # QC the unmapped BAM 
    call CollectQualityYieldMetrics {
      input:
        PICARD=picard,
        in_bam = unmapped_bam,
        metrics_filename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmapped.quality_yield_metrics"
    }

    # Map reads to reference
    call SamToFastqAndBwaMem {
      input:
        SAMTOOLS=samtools,
        PICARD=picard,
        in_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        out_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa
     }

    # Merge original uBAM and BWA-aligned BAM 
    call MergeBamAlignment {
      input:
        PICARD=picard,
        unmapped_bam = unmapped_bam,
        aligned_bam = SamToFastqAndBwaMem.out_bam,
        out_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }

    # QC the aligned but unsorted readgroup BAM
    # No reference needed as the input here is unsorted; providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment.out_bam,
        out_bam_prefix = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".readgroup"
    }

    # Sort and fix tags in the merged BAM
    call SortAndFixTags as SortAndFixReadGroupBam {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment.out_bam,
        out_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".sorted",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
    
    # Validate the aligned and sorted readgroup BAM
    # This is called to help in finding problems early.
    # If considered too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile {
      input:
        PICARD=picard,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        in_bam = SortAndFixReadGroupBam.out_bam,
        in_bam_index = SortAndFixReadGroupBam.out_bam_index,
        report_filename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".validation_report"
    }

  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      PICARD=picard,
      in_bams = MergeBamAlignment.out_bam,
      out_bam_basename = base_file_name_normal + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name_normal + ".duplicate_metrics"
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam {
    input:
      PICARD=picard,
      in_bam = MarkDuplicates.out_bam,
      out_bam_basename = base_file_name_normal + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  # Check identity of fingerprints across readgroups
  call CrossCheckFingerprints {
    input:
      PICARD=picard,
      in_bams = SortAndFixSampleBam.out_bam,
      in_bam_indexes = SortAndFixSampleBam.out_bam_index,
      haplotype_database_file = haplotype_database_file,
      metrics_filename = base_file_name_normal + ".crosscheck"
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }
  
  # Estimate level of cross-sample contamination
  call CheckContamination {
    input:
      verifyBamID=verifyBamID,
      in_bam = SortAndFixSampleBam.out_bam,
      in_bam_index = SortAndFixSampleBam.out_bam_index,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      output_prefix = base_file_name_normal + ".preBqsr"
  }
  
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        GATK=gatk,
        in_bam = SortAndFixSampleBam.out_bam,
        in_bam_index = SortAndFixSampleBam.out_bam_index,
        recalibration_report_filename = base_file_name_normal + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }  
  }  
  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports {
    input:
      GATK=gatk,
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name_normal + ".recal_data.csv"
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        GATK4=gatk4,
        in_bam = SortAndFixSampleBam.out_bam,
        in_bam_index = SortAndFixSampleBam.out_bam_index,
        out_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      PICARD=picard,
      in_bams = ApplyBQSR.recalibrated_bam,
      out_bam_basename = base_file_name_normal
  }
  
  # QC the final BAM (consolidated after scattered BQSR)
  call CollectReadgroupBamQualityMetrics {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      out_bam_prefix = base_file_name_normal + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  # Validate the final BAM 
  call ValidateSamFile as ValidateAggregatedSamFile {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      report_filename = base_file_name_normal + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  
  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      out_bam_prefix = base_file_name_normal,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  
  # Check the sample BAM fingerprint against the sample array 
  call CheckFingerprint {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      haplotype_database_file = haplotype_database_file,
      genotypes = fingerprint_genotypes_file,
      output_basename = base_file_name_normal,
      sample = sample_name + '_normal'
  }
  
  # QC the sample WGS metrics (stringent thresholds)
  call CollectWgsMetrics {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      metrics_filename = base_file_name_normal + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list
  }
  
  # QC the sample raw WGS metrics (common thresholds)
  call CollectRawWgsMetrics {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      metrics_filename = base_file_name_normal + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list
  }
  
  # Generate a checksum per readgroup in the final BAM
  call CalculateReadGroupChecksum {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5"
  }
  
  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram {
    input:
      seq_cache_populate=seq_cache_populate,
      SAMTOOLS=samtools,
      in_bam = GatherBamFiles.out_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = base_file_name_normal
  }

  # Convert the CRAM back to BAM to check that the conversions do not introduce errors
  call CramToBam {
    input:
      SAMTOOLS=samtools,
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_fasta_index = ref_fasta_index,
      cram_file = ConvertToCram.output_cram,
      output_basename = base_file_name_normal + ".roundtrip"
  }

  # Validate the roundtripped BAM
  call ValidateSamFile as ValidateBamFromCram {
    input:
      PICARD=picard,
      in_bam = CramToBam.out_bam,
      in_bam_index = CramToBam.out_bam_index,
      report_filename = base_file_name_normal + ".bam.roundtrip.validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      max_output = 1000000000,
      ignore = ["null"]
  }
  
  # Call variants in parallel over WGS calling intervals
  scatter (subInterval in scattered_calling_intervals) {
  
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        GATK=gatk,
        contamination = CheckContamination.contamination,
        in_bam = GatherBamFiles.out_bam,
        in_bam_index = GatherBamFiles.out_bam_index,
        interval_list = subInterval,
        gvcf_basename = base_file_name_normal,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
     }
  }
  
  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs {
    input:
      PICARD=picard,
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = final_gvcf_name_normal
  }
  
  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF {
    input:
      GATK=gatk,
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list
  }
  
  # QC the GVCF
  call CollectGvcfCallingMetrics {
    input:
      PICARD=picard,
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = base_file_name_normal,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list
  }























  String final_gvcf_name_tumor = base_file_name_tumor + final_gvcf_ext


  ###### Here the workflow is adapted to .fastq files:
  # Decompress and split the files into chunks, return an array of .fastq files:
  call UnzipAndSplit as UnzipAndSplit_tumor {
      input:
        PIGZ=pigz,
        input_fastqR1 = rawdata_normal_fastqR1,
        input_fastqR2 = rawdata_normal_fastqR2,
        sample_name = sample_name + '_tumor'
  }


  Array[Pair[File, File]] fastqR1R2_chunks_tumor = zip(UnzipAndSplit_tumor.R1_splits, UnzipAndSplit_tumor.R2_splits)
  scatter (fastq_chunk_tumor in fastqR1R2_chunks_tumor) {

    call TrimReads as TrimReads_tumor {
        input:
          TRIMMOMATIC=trimmomatic,
          input_fastqR1 = fastq_chunk_tumor.left,
          input_fastqR2 = fastq_chunk_tumor.right,
          basenameR1 = basename(fastq_chunk_tumor.left, ".fastq"),
          basenameR2 = basename(fastq_chunk_tumor.right, ".fastq")
    }

    call FastqToBam as FastqToBam_tumor {
        input:
          PICARD=picard,
          sample_name = sample_name + '_tumor',
          input_fastqR1 = TrimReads_tumor.output_R1,
          input_fastqR2 = TrimReads_tumor.output_R2
    }

#  }

  ######


  # Align flowcell-level unmapped input bams in parallel
#  scatter (unmapped_bam in flowcell_unmapped_bams) {
#  scatter (unmapped_bam in FastqToBam.out_bam) {  
    # Because of a wdl/cromwell bug this is not currently valid so we have to sub(sub()) in each task
    # String base_name = sub(sub(unmapped_bam, "gs://.*/", ""), unmapped_bam_suffix + "$", "")

    String sub_strip_path = "/home/projects/dp_00005/data/cromwell_test/.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"

    # QC the unmapped BAM 
    call CollectQualityYieldMetrics as CollectQualityYieldMetrics_tumor {
      input:
        PICARD=picard,
        in_bam = FastqToBam_tumor.out_bam,
        metrics_filename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmapped.quality_yield_metrics"
    }

    call STAR_Map {
      input: STAR=star,
        STARindexDir=premade_STARindexDir,
        sample_name = sample_name + '_tumor',
        input_fastqR1 = TrimReads_tumor.output_R1,
        input_fastqR2 = TrimReads_tumor.output_R2
    }

#    # Map reads to reference
#    call SamToFastqAndBwaMem {
#      input:
#        SAMTOOLS=samtools,
#        PICARD=picard,
#        in_bam = FastqToBam_tumor.out_bam,
#        out_bam_basename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmerged",
#        ref_fasta = ref_fasta,
#        ref_fasta_index = ref_fasta_index,
#        ref_dict = ref_dict,
#        ref_alt = ref_alt,
#        ref_bwt = ref_bwt,
#        ref_amb = ref_amb,
#        ref_ann = ref_ann,
#        ref_pac = ref_pac,
#        ref_sa = ref_sa
#     }


    # Merge original uBAM and BWA-aligned BAM 
    call MergeBamAlignment as MergeBamAlignment_tumor {
      input:
        PICARD=picard,
        unmapped_bam = FastqToBam_tumor.out_bam,
        aligned_bam = STAR_Map.out_bam,
        out_bam_basename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }

    # QC the aligned but unsorted readgroup BAM
    # No reference needed as the input here is unsorted; providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics_tumor {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment_tumor.out_bam,
        out_bam_prefix = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".readgroup"
    }

    # Sort and fix tags in the merged BAM
    call SortAndFixTags as SortAndFixReadGroupBam_tumor {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment_tumor.out_bam,
        out_bam_basename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".sorted",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
    
    # Validate the aligned and sorted readgroup BAM
    # This is called to help in finding problems early.
    # If considered too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile_tumor {
      input:
        PICARD=picard,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        in_bam = SortAndFixReadGroupBam_tumor.out_bam,
        in_bam_index = SortAndFixReadGroupBam_tumor.out_bam_index,
        report_filename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".validation_report"
    }

  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates as MarkDuplicates_tumor {
    input:
      PICARD=picard,
      in_bams = MergeBamAlignment_tumor.out_bam,
      out_bam_basename = base_file_name_tumor + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name_tumor + ".duplicate_metrics"
  }

  # Use SplitNCigarReads for best practices on RNAseq data:
  call SplitNCigarReads {
    input: GATK=gatk,
      sample_name = sample_name + '_tumor',
      ref_fasta=ref_fasta,
      in_bam=MarkDuplicates_tumor.out_bam
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam_tumor {
    input:
      PICARD=picard,
      in_bam = SplitNCigarReads.out_bam,
      out_bam_basename = base_file_name_tumor + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  # Check identity of fingerprints across readgroups
  call CrossCheckFingerprints as CrossCheckFingerprints_tumor {
    input:
      PICARD=picard,
      in_bams = SortAndFixSampleBam_tumor.out_bam,
      in_bam_indexes = SortAndFixSampleBam_tumor.out_bam_index,
      haplotype_database_file = haplotype_database_file,
      metrics_filename = base_file_name_tumor + ".crosscheck"
  }

  # Create list of sequences for scatter-gather parallelization 
  call CreateSequenceGroupingTSV as CreateSequenceGroupingTSV_tumor {
    input:
      ref_dict = ref_dict
  }
  
  # Estimate level of cross-sample contamination
  call CheckContamination as CheckContamination_tumor {
    input:
      verifyBamID=verifyBamID,
      in_bam = SortAndFixSampleBam_tumor.out_bam,
      in_bam_index = SortAndFixSampleBam_tumor.out_bam_index,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      output_prefix = base_file_name_tumor + ".preBqsr"
  }
  
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV_tumor.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator as BaseRecalibrator_tumor {
      input:
        GATK=gatk,
        in_bam = SortAndFixSampleBam_tumor.out_bam,
        in_bam_index = SortAndFixSampleBam_tumor.out_bam_index,
        recalibration_report_filename = base_file_name_tumor + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }  
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports as GatherBqsrReports_tumor {
    input:
      GATK=gatk,
      input_bqsr_reports = BaseRecalibrator_tumor.recalibration_report,
      output_report_filename = base_file_name_tumor + ".recal_data.csv"
  }

  scatter (subgroup in CreateSequenceGroupingTSV_tumor.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR as ApplyBQSR_tumor {
      input:
        GATK4=gatk4,
        in_bam = SortAndFixSampleBam_tumor.out_bam,
        in_bam_index = SortAndFixSampleBam_tumor.out_bam_index,
        out_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports_tumor.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      PICARD=picard,
      in_bams = ApplyBQSR_tumor.recalibrated_bam,
      out_bam_basename = base_file_name_tumor
  }
  
  # QC the final BAM (consolidated after scattered BQSR)
  call CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      out_bam_prefix = base_file_name_tumor + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }

  # Validate the final BAM 
  call ValidateSamFile as ValidateAggregatedSamFile_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      report_filename = base_file_name_tumor + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  
  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics as CollectAggregationMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      out_bam_prefix = base_file_name_tumor,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index
  }
  
  # Check the sample BAM fingerprint against the sample array 
  call CheckFingerprint as CheckFingerprint_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      haplotype_database_file = haplotype_database_file,
      genotypes = fingerprint_genotypes_file,
      output_basename = base_file_name_tumor,
      sample = sample_name + '_tumor'
  }
  
  # QC the sample WGS metrics (stringent thresholds)
  call CollectWgsMetrics as CollectWgsMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      metrics_filename = base_file_name_tumor + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list
  }
  
  # QC the sample raw WGS metrics (common thresholds)
  call CollectRawWgsMetrics as CollectRawWgsMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      metrics_filename = base_file_name_tumor + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list
  }
  
  # Generate a checksum per readgroup in the final BAM
  call CalculateReadGroupChecksum as CalculateReadGroupChecksum_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles.out_bam,
      in_bam_index = GatherBamFiles.out_bam_index,
      read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5"
  }
  
  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram as ConvertToCram_tumor {
    input:
      seq_cache_populate=seq_cache_populate,
      SAMTOOLS=samtools,
      in_bam = GatherBamFiles.out_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = base_file_name_tumor
  }

  # Convert the CRAM back to BAM to check that the conversions do not introduce errors
  call CramToBam as CramToBam_tumor {
    input:
      SAMTOOLS=samtools,
      ref_fasta = ref_fasta,
      ref_dict = ref_dict,
      ref_fasta_index = ref_fasta_index,
      cram_file = ConvertToCram_tumor.output_cram,
      output_basename = base_file_name_tumor + ".roundtrip"
  }

  # Validate the roundtripped BAM
  call ValidateSamFile as ValidateBamFromCram_tumor {
    input:
      PICARD=picard,
      in_bam = CramToBam_tumor.out_bam,
      in_bam_index = CramToBam_tumor.out_bam_index,
      report_filename = base_file_name_tumor + ".bam.roundtrip.validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      max_output = 1000000000,
      ignore = ["null"]
  }
  
  # Call variants in parallel over WGS calling intervals
  scatter (subInterval in scattered_calling_intervals) {
  
    # Generate GVCF by interval
    call HaplotypeCaller_RNA as HaplotypeCaller_tumor {
      input:
        GATK=gatk,
        contamination = CheckContamination_tumor.contamination,
        in_bam = GatherBamFiles.out_bam,
        in_bam_index = GatherBamFiles.out_bam_index,
        interval_list = subInterval,
        gvcf_basename = base_file_name_tumor,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index
     }
  }
  
  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs as MergeVCFs_tumor {
    input:
      PICARD=picard,
      input_vcfs = HaplotypeCaller_tumor.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller_tumor.output_gvcf_index,
      output_vcf_name = final_gvcf_name_tumor
  }
  
  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF as ValidateGVCF_tumor {
    input:
      GATK=gatk,
      input_vcf = MergeVCFs_tumor.output_vcf,
      input_vcf_index = MergeVCFs_tumor.output_vcf_index,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list
  }
  
  # QC the GVCF
  call CollectGvcfCallingMetrics as CollectGvcfCallingMetrics_tumor {
    input:
      PICARD=picard,
      input_vcf = MergeVCFs_tumor.output_vcf,
      input_vcf_index = MergeVCFs_tumor.output_vcf_index,
      metrics_basename = base_file_name_tumor,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list
  }





















  # Outputs that will be retained when execution is complete  
  output {
    CollectQualityYieldMetrics.*
    CollectQualityYieldMetrics_tumor.*
    ValidateReadGroupSamFile.*
    ValidateReadGroupSamFile_tumor.*
    CollectReadgroupBamQualityMetrics.*
    CollectReadgroupBamQualityMetrics_tumor.*
    CollectUnsortedReadgroupBamQualityMetrics.*
    CollectUnsortedReadgroupBamQualityMetrics_tumor.*
    CrossCheckFingerprints.*
    CrossCheckFingerprints_tumor.*
    ValidateBamFromCram.*
    ValidateBamFromCram_tumor.*
    CalculateReadGroupChecksum.*
    CalculateReadGroupChecksum_tumor.*
    ValidateAggregatedSamFile.*
    ValidateAggregatedSamFile_tumor.*
    CollectAggregationMetrics.*
    CollectAggregationMetrics_tumor.*
    CheckFingerprint.*
    CheckFingerprint_tumor.*
    CollectWgsMetrics.*
    CollectWgsMetrics_tumor.*
    CollectRawWgsMetrics.*
    CollectRawWgsMetrics_tumor.*
    CheckContamination.*
    CheckContamination_tumor.*
    CollectGvcfCallingMetrics.*
    CollectGvcfCallingMetrics_tumor.*
    MarkDuplicates.duplicate_metrics
    MarkDuplicates_tumor.duplicate_metrics
    GatherBqsrReports.*
    GatherBqsrReports_tumor.*
    ConvertToCram.*
    ConvertToCram_tumor.*
    MergeVCFs.*
    MergeVCFs_tumor.*
    } 
}
