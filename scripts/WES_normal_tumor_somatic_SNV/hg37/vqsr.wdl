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
task MuTect2 {
  String sample_name
  File in_bam_tumor
  File in_bai_tumor
  File in_bam_normal
  File in_bai_normal
  String sample_name_tumor
  String sample_name_normal
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  File transcript_intervals
  Float? contamination
  Int cpu=28
  File GATK4_LAUNCH
  File dbSNP_vcf
  File dbSNP_vcf_index

  command {
    ${GATK4_LAUNCH} --javaOptions "-Xmx10g" Mutect2 \
     -R ${ref_fasta} \
     -I ${in_bam_tumor} \
     -tumor ${sample_name_tumor} \
     -I ${in_bam_normal} \
     -normal ${sample_name_tumor} \
     --dbsnp ${dbSNP_vcf} \
     --dontUseSoftClippedBases \
     -L ${transcript_intervals} \
     -O ${sample_name}_MuTect2.vcf.gz \
     --contamination_fraction_to_filter ${default=0 contamination}
  }
  runtime {
    cpu: cpu
  }
  output {
    File output_vcf = "${sample_name}_MuTect2.vcf.gz"
    File output_vcf_index = "${sample_name}_MuTect2.vcf.gz.tbi"
  }
}

task STAR_Map {
  File STAR
  String STARindexDir
  String sample_name
  File input_fastqR1
  File input_fastqR2
#  String suffix="Aligned.sortedByCoord.out"
  # This is STAR's suffix and cannot change:
  String suffix="Aligned.out"
  Int cpu=28

  command {
    ln -s ${STARindexDir} GenomeDir
    # Use default mode --genomeDir=./GenomeDir/
    ${STAR} --readFilesIn ${input_fastqR1} ${input_fastqR2} \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM Unsorted \
        --outFileNamePrefix ${sample_name} \
        --twopassMode Basic \
        --outSAMmapqUnique 60 \
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
  File in_bai
  String suffix="_splitNcigar"
  Int cpu=2

  command {
    # Run options are set to follow the GATK best practice for RNAseq. data.
    /home/projects/dp_00005/apps/src/gatk-4.beta.5/gatk-launch \
      SplitNCigarReads \
      -R ${ref_fasta} \
      -I ${in_bam} \
      -O ${sample_name}${suffix}.bam
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
  File R
  File Rscript
  File samtools 
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
    File out_bai = "${out_bam_basename}.bai"
    File out_bam_md5 = "${out_bam_basename}.bam.md5"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  File in_bam
  String out_bam_prefix
  Int cpu=1
  File PICARD
  File R
  File Rscript 
  File samtools 

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
  File in_bai
  String out_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD
  File R
  File Rscript
  File samtools 

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
  File in_bai
  String out_bam_prefix
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD
  File R 
  File Rscript
  File samtools 

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
    File insert_size_histogram_pdf = "${out_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "${out_bam_prefix}.insert_size_metrics"
  }
}

# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  Array[File] in_bams
  Array[File] in_baies
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
  File in_bai
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
  File PYTHON2
  File ref_dict
  Int cpu=1

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    ${PYTHON2} <<CODE
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
  File in_bai
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
  String? U_option

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
      -L $rand.intervals \
      ${U_option}
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
  File in_bai
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
    File out_bai = "${out_bam_basename}.bai"
    File out_bam_md5 = "${out_bam_basename}.bam.md5"
  }
}

# Validate the output bam file
task ValidateSamFile {
  File in_bam
  File in_bai
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
  File in_bai
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD
  File R 
  File Rscript 
  File samtools 

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
  File in_bai
  String metrics_filename
  File wgs_coverage_interval_list
  File ref_fasta
  File ref_fasta_index
  Int cpu=1
  File PICARD
  File R 
  File Rscript
  File samtools

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
  File in_bai
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
  File in_bai
  File contamination_sites_vcf
  File contamination_sites_vcf_index
  String output_prefix
  Int cpu=1
  File verifyBamID
  File PYTHON3

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

    ${PYTHON3} <<CODE
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
  File in_bai
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int cpu=28
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
      -nct ${cpu} \
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
  File in_bai
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  Int cpu=28
  File GATK

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx20000m \
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
      -nct 5 \
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

# Convert combined gVCF from MergeVCFs and HaplotypeCaller 
task GenotypeGVCF {
  File GATK
  File ref_fasta
  File ref_fasta_index
  File ref_dict 
  File input_gvcf 
  File input_gvcf_index
  String output_basename
  Int cpu=1
  
  command {
    java -jar ${GATK} \
     -T GenotypeGVCFs \
     -R ${ref_fasta} \
     --variant ${input_gvcf} \
     -o ${output_basename}.vcf
  }
  output {
     File output_vcf = "${output_basename}.vcf"
     File output_vcf_index = "${output_basename}.vcf.idx"
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
  File R 
  File Rscript
  File samtools

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
    File out_bai = "${output_basename}.bai"
  }
}

# FILTERING TASKS # 

task VariantRecalibrator_SNP {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String output_basename
  Int cpu=28
  File hapmap
  File omni
  File Genomes
  File dbSNP_vcf
  File GATK
  File input_vcf
  File input_vcf_index

  command {
    java -jar ${GATK} \
       -T VariantRecalibrator \
       -R ${ref_fasta} \
       -input ${input_vcf} \
       -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
       -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
       -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${Genomes} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP_vcf} \
       -an QD \
       -an FS \
       -an SOR \
       -an MQ \
       -an MQRankSum \
       -an ReadPosRankSum \
       -mode SNP \
       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
       -recalFile SNP.${output_basename}.recal \
       -tranchesFile SNP.${output_basename}.tranches
  }
  runtime {
    cpu: cpu
  }
  output {
    File recalFile = "SNP.${output_basename}.recal"
    File recalFile_index = "SNP.${output_basename}.recal.idx"
    File tranchesFile = "SNP.${output_basename}.tranches"
  }
}

## In lack of programming skills for handling resource input, I define the task for INDEL VariantRecalibration again - in the future it would be cleaner to use the same task twice
task VariantRecalibrator_INDEL {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String output_basename
  Int cpu=28
  File mills 
  File dbSNP_vcf
  File GATK
  File input_vcf
  File input_vcf_index

  command {
    java -jar ${GATK} \
       -T VariantRecalibrator \
       -R ${ref_fasta} \
       -input ${input_vcf} \
       -resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills} \
       -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP_vcf} \
       -an QD \
       -an FS \
       -an SOR \
       -an MQ \
       -an MQRankSum \
       -an ReadPosRankSum \
       -mode INDEL \
       --maxGaussians 4 \
       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
       -recalFile INDEL.${output_basename}.recal \
       -tranchesFile INDEL.${output_basename}.tranches
  }
  runtime {
    cpu: cpu
  }
  output {
    File recalFile = "INDEL.${output_basename}.recal"
    File recalFile_index = "INDEL.${output_basename}.recal.idx"
    File tranchesFile = "INDEL.${output_basename}.tranches"
  }
}

# Calibrate variant calls
task ApplyVQSR {
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String output_basename
  Int cpu=1
  File GATK
  File input_vcf
  File input_vcf_index
  File recalFile
  File recalFile_index
  File tranchesFile
  String MODE
  Float tranchesLevel

  command {
    java -jar ${GATK} \
      -T ApplyRecalibration \
      -R ${ref_fasta} \
      -input ${input_vcf} \
      -mode ${MODE} \
      --ts_filter_level ${tranchesLevel} \
      -recalFile ${recalFile} \
      -tranchesFile ${tranchesFile} \
      -o ${MODE}_recalibrated.${output_basename}.vcf
  }
  runtime {
    cpu: cpu
  }
  output {
    File output_vcf = "${MODE}_recalibrated.${output_basename}.vcf"
    File output_vcf_index = "${MODE}_recalibrated.${output_basename}.vcf.idx"
  }
}

# filtering RNA Seq variants 
task VariantFiltration {
  File ref_fasta 
  File ref_fasta_index
  File ref_dict
  String output_basename
  Int cpu=1
  File GATK
  File input_vcf
  File input_vcf_index
  Int window_size
  Float QD
  Float FS

  command {
    java -jar ${GATK} \
       -T VariantFiltration \
       -R ${ref_fasta} \
       -V ${input_vcf} \
       -window ${window_size} \
       -cluster 3 \
       -filterName FS \
       -filter "FS > ${FS}" \
       -filterName QD \
       -filter "QD < ${QD}" \
       -o ${output_basename}_filtered.vcf
  }
  runtime {
    cpu: cpu
  }
  output {
    File output_vcf = "${output_basename}_filtered.vcf"
    File output_vcf_index = "${output_basename}_filtered.vcf.idx"
  }
}   

# WORKFLOW DEFINITION
workflow WGS_normal_RNAseq_tumor_SNV_wf {

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
  String premade_STARindexDir

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File transcript_intervals

  # Filtering resources: 
  File hapmap
  File mills
  File Genomes
  File omni

  String recalibrated_bam_basename = base_file_name_normal + ".aligned.duplicates_marked.recalibrated"

  # Tools:
  File picard
  File gatk
  File gatk4
  File gatk4_launch
  File python2
  File python3
  File samtools
  File bwa
  File verifyBamID
  File seq_cache_populate
  File pigz
  File trimmomatic
  File star
  File gcc
  File R
  File Rscript
  File OracleJava 

  String bwa_commandline = bwa + " mem -K 100000000 -p -v 3 -t 28 -Y $bash_ref_fasta"

 
  String final_gvcf_name_normal = base_file_name_normal + final_gvcf_ext
   
  # Filtering input 
  File input_vcf
  File input_vcf_index

  # following two lines are used for both DNA and RNA part 
  String sub_strip_path = "/home/projects/dp_00005/data/cromwell_test/.*/"
  String sub_strip_unmapped = unmapped_bam_suffix + "$"

#####################
#### DNA only part ###
######################

 # Filter HaplotypeCaller output using VQSR 
 # Step 1: SNP recalibration 
 call VariantRecalibrator_SNP as VariantRecalibrator_SNP_normal {
   input: 
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict, 
    output_basename = base_file_name_normal, 
    hapmap=hapmap, 
    omni=omni,
    Genomes=Genomes,
    dbSNP_vcf=dbSNP_vcf,
    GATK=gatk, 
    input_vcf=input_vcf,
    input_vcf_index=input_vcf_index
  }

 # Step 2: Apply SNP recalibration 
 call ApplyVQSR as ApplyVQSR_SNP_normal {
   input:
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict, 
    output_basename = base_file_name_normal,
    GATK=gatk,
    input_vcf=input_vcf,
    input_vcf_index=input_vcf_index,
    recalFile=VariantRecalibrator_SNP_normal.recalFile,
    recalFile_index=VariantRecalibrator_SNP_normal.recalFile_index,
    tranchesFile=VariantRecalibrator_SNP_normal.tranchesFile,
    tranchesLevel=99.0,
    MODE="SNP"
  }
 
 # Step 3: Indel recalibration
 call VariantRecalibrator_INDEL as VariantRecalibrator_INDEL_normal {
   input:
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict, 
    output_basename = base_file_name_normal,
    mills=mills,
    dbSNP_vcf=dbSNP_vcf,
    GATK=gatk,
    input_vcf=ApplyVQSR_SNP_normal.output_vcf,
    input_vcf_index=ApplyVQSR_SNP_normal.output_vcf_index
  }
  
 # Step 4: Apply Indel recalibration 
  call ApplyVQSR as ApplyVQSR_INDEL_normal {
   input:
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index,
    ref_dict = ref_dict,
    output_basename = base_file_name_normal,
    GATK=gatk,
    input_vcf=ApplyVQSR_SNP_normal.output_vcf,
    input_vcf_index=ApplyVQSR_SNP_normal.output_vcf_index,
    recalFile=VariantRecalibrator_INDEL_normal.recalFile,
    recalFile_index=VariantRecalibrator_INDEL_normal.recalFile_index,
    tranchesFile=VariantRecalibrator_INDEL_normal.tranchesFile,
    tranchesLevel=99.0,
    MODE = "INDEL"
  }

#########################
### DNA only part end ###
#########################


  # Outputs that will be retained when execution is complete  
  output {
    VariantRecalibrator_SNP_normal.*
    VariantRecalibrator_INDEL_normal.*
    ApplyVQSR_SNP_normal.*
    ApplyVQSR_INDEL_normal.*
    } 
}
