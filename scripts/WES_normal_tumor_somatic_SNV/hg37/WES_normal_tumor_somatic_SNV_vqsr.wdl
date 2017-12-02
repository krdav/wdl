## Based on Broad Institute germline SNP and indel discovery workflow:
## https://github.com/broadinstitute/wdl/blob/develop/scripts/broad_pipelines/PublicPairedSingleSampleWf_170412.wdl
## Copyright Broad Institute, 2017
##
## The workflow has been adapted towards calling somatic tumor variants using a pair
## of normal/tumor sample from either WES or WGS. The WDL pipeline follows GATK's
## Best Practices (as of June 2016).
##
## Requirements/expectations:
## - Human WES/WGS paired-end sequencing data in fastq format for both a normal
##   and a tumor sample.
##
## Runtime parameters are optimized for DTU's Computerome Google HPC platform.
##
## LICENSING:
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script.


########################
### TASK DEFINITIONS ###
########################

# Build VQSR model:
task BuildVQSRModel {
  File ref_dict
  File ref_fa 
  File ref_idx
  Int cpu=1
  File GATK
  File in_vcf
  File in_vcf_idx
  File wgs_calling_interval_list
  String out_basename
  String mode
  Array[String] annotations
  Array[Float] tranches
  Array[String] resources
  Array[File] resource_files
  Array[File] resource_indices

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar ${GATK} \
      -T VariantRecalibrator \
      -R ${ref_fa} \
      -input ${in_vcf} \
      -L ${interval_list} \
      -resource:${sep=' -resource:' resources} \
      -an ${sep=' -an ' annotations} \
      -mode ${mode} \
      -tranche ${sep=' -tranche ' tranches} \
      -recalFile ${out_basename}.${mode}.recal \
      -tranchesFile ${out_basename}.${mode}.tranches \
      -rscriptFile ${out_basename}.${mode}.plots.R
  }

  runtime {
    cpu: cpu
  }

  output {
    File recal_file = "${out_basename}.${mode}.recal"
    File recal_file_idx = "${out_basename}.${mode}.recal.idx"
    File tranches_file = "${out_basename}.${mode}.tranches"
    File rscript_file = "${out_basename}.${mode}.plots.R"
  }
}




# Apply recalibration:
task ApplyRecalibrationFilter {
  File ref_dict
  File ref_fa 
  File ref_idx
  Int cpu=1
  File GATK
  File in_vcf
  File in_vcf_idx
  String out_basename
  File wgs_calling_interval_list
  File recal_file
  File recal_file_idx
  String mode
  File tranches_file
  Float filter_level

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar ${GATK} \
      -T ApplyRecalibration \
      -R ${ref_fa} \
      -input ${in_vcf} \
      -L ${interval_list} \
      -mode ${mode} \
      --ts_filter_level ${filter_level} \
      -recalFile ${recal_file} \
      -tranchesFile ${tranches_file} \
      -o ${out_vcf_name}
  }

  runtime {
    cpu: cpu
  }

  output {
    File out_vcf = "${out_basename}.g.vcf.gz"
    File out_vcf_idx = "${out_basename}.g.vcf.gz.tbi"
  }
}































# Sequence-context dependent artifacts filtering of MuTect2 variants:
task FilterByOrientationBias {
  String sample_name
  File in_vcf
  File in_vcf_idx
  Int cpu=1
  File GATK4_LAUNCH
  File pre_adapter_detail_metrics_tumor
  String pre_adapter_detail_metrics_tumor_base = basename(pre_adapter_detail_metrics_tumor)

  command {
    # Nasty sed replacing hack because of: https://github.com/broadinstitute/gatk/issues/3030
    sed -r "s/picard\.analysis\.artifacts\.SequencingArtifactMetrics\\\$PreAdapterDetailMetrics/org\.broadinstitute\.hellbender\.tools\.picard\.analysis\.artifacts\.SequencingArtifactMetrics\$PreAdapterDetailMetrics/g" "${pre_adapter_detail_metrics_tumor}" > "gatk_${pre_adapter_detail_metrics_tumor_base}"

    ${GATK4_LAUNCH} --javaOptions "-Xmx10g" FilterByOrientationBias \
      -V ${in_vcf} \
      -O ${sample_name}_MuTect2_filtered2.vcf.gz \
      -P "gatk_${pre_adapter_detail_metrics_tumor_base}" \
      -A 'G/T' -A 'C/T'
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_vcf = "${sample_name}_MuTect2_filtered2.vcf.gz"
    File out_vcf_idx = "${sample_name}_MuTect2_filtered2.vcf.gz.tbi"
  }
}

# Simple filtering of MuTect2 variants:
task FilterMutectCalls {
  String sample_name
  File in_vcf
  File in_vcf_idx
  Float? contamination
  Int cpu=1
  File GATK4_LAUNCH
  File dbSNP_vcf
  File dbSNP_vcf_idx

  command {
    ${GATK4_LAUNCH} --javaOptions "-Xmx10g" FilterMutectCalls \
      -V ${in_vcf} \
      -O ${sample_name}_MuTect2_filtered.vcf.gz \
      --dbsnp ${dbSNP_vcf} \
      --contamination_fraction_to_filter ${default=0 contamination} \
      --normal_artifact_lod 4.0 \
      --tumor_lod 4.0 \
      --base_quality_score_threshold 20 \
      --min_base_quality_score 20 \
      --max_germline_posterior 0.001 \
      --pcr_indel_model HOSTILE \
      --standard_min_confidence_threshold_for_calling 20 \
      --strandArtifactAlleleFraction 0.1 \
      --strandArtifactPosteriorProbability 0.1 \
      --uniqueAltReadCount 5
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_vcf = "${sample_name}_MuTect2_filtered.vcf.gz"
    File out_vcf_idx = "${sample_name}_MuTect2_filtered.vcf.gz.tbi"
  }
}

# MuTect2 somatic variant calling with normal/tumor sample input:
task MuTect2 {
  String sample_name
  File in_bam_tumor
  File in_bai_tumor
  File in_bam_normal
  File in_bai_normal
  String sample_name_tumor
  String sample_name_normal
  File ref_dict
  File ref_fa
  File ref_idx
  File transcript_intervals
  Float? contamination
  Int cpu=28
  File GATK4_LAUNCH
  File gnomad_exome_vcf
  File gnomad_exome_vcf_idx
  File dbSNP_vcf
  File dbSNP_vcf_idx

  command {
    ${GATK4_LAUNCH} --javaOptions "-Xms5g -Xmx100g" Mutect2 \
      -R ${ref_fa} \
      -I ${in_bam_tumor} \
      -tumor ${sample_name_tumor} \
      -I ${in_bam_normal} \
      -normal ${sample_name_tumor} \
      --dbsnp ${dbSNP_vcf} \
      --dontUseSoftClippedBases \
      -L ${transcript_intervals} \
      -O ${sample_name}_MuTect2.vcf.gz \
      --germline_resource ${gnomad_exome_vcf} \
      --contamination_fraction_to_filter ${default=0 contamination}
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_vcf = "${sample_name}_MuTect2.vcf.gz"
    File out_vcf_idx = "${sample_name}_MuTect2.vcf.gz.tbi"
  }
}

# First split the fastq files into chunks of 10 million reads
task UnzipAndSplit {
  File in_fastqR1
  File in_fastqR2
  File PIGZ
  String sample_name
  Int cpu=2

  # Uncompress while splitting:
  command {
    ${PIGZ} -dc -p 2 ${in_fastqR1} | split -l 40000000 --additional-suffix=".fastq" - "${sample_name}_1_" &
    ${PIGZ} -dc -p 2 ${in_fastqR2} | split -l 40000000 --additional-suffix=".fastq" - "${sample_name}_2_" &
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

# Read trimming (can be omitted because soft clipping by the aligner should deal with this)
task TrimReads {
  File in_fastqR1
  File in_fastqR2
  String basenameR1
  String basenameR2
  File TRIMMOMATIC
  File adapters = '/services/tools/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa'
  Int cpu=4  # More CPUs will not decrease runtime substantially

  command {
    java -Xmx80g \
      -jar ${TRIMMOMATIC} \
      PE \
      -threads ${cpu} \
      -phred33 \
      ${in_fastqR1} ${in_fastqR2} ${basenameR1}_trim.fastq ${basenameR1}_trim_unpaired.fastq ${basenameR2}_trim.fastq ${basenameR2}_trim_unpaired.fastq \
      ILLUMINACLIP:${adapters}:2:30:10
  }
  runtime {
    cpu: cpu
  }
  output {
    File out_R1 = "${basenameR1}_trim.fastq"
    File out_R2 = "${basenameR2}_trim.fastq"
  }
}

# Convert two paired-end fastq files into a single bam file
task FastqToBam {
  File in_fastqR1
  File in_fastqR2
  String basenameR1 = basename(in_fastqR1, ".fastq")
  String basename = sub(basenameR1, '_1', '')
  File PICARD
  String sample_name
  Int cpu=2

  command {
    java -Xmx8G \
      -jar ${PICARD} FastqToSam \
      FASTQ=${in_fastqR1} \
      FASTQ2=${in_fastqR2} \
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

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
task SamToFastqAndBwaMem {
  File in_bam
  String bwa_commandline
  String out_bam_basename
  File ref_fa
  File ref_idx
  File ref_dict
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
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
    bash_ref_fa=${ref_fa}
    # if ref_amb has data in it, we proceed
    if [ -s ${ref_amb} ]; then
      java -Xmx3000m \
        -jar ${PICARD} \
        SamToFastq \
        INPUT=${in_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ${bwa_commandline} /dev/stdin -  2> >(tee ${out_bam_basename}.bwa.stderr.log >&2) | \
      ${SAMTOOLS} view -1 - > ${out_bam_basename}.bam

      grep -m1 "read .* ALT contigs" ${out_bam_basename}.bwa.stderr.log #| \
      # grep -v "read 0 ALT contigs"

    # else ref_amb is empty or could not be found, so we bail out
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
  File ref_fa
  File ref_idx
  File ref_dict
  Int cpu=1
  File PICARD

  command {
    # set the bash variable needed for the command-line
    bash_ref_fa=${ref_fa}
    java -Xmx2500m \
      -jar ${PICARD} \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${out_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fa} \
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
  File ref_fa
  File ref_idx
  Int cpu=2
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
      REFERENCE_SEQUENCE=${ref_fa}
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
  File ref_fa
  File ref_idx
  Int cpu=1
  File PICARD

  command {
    java -Xmx5000m \
      -jar ${PICARD} \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      REFERENCE_SEQUENCE=${ref_fa} \
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
  File ref_fa
  File ref_idx
  Int cpu=1
  File PICARD

  command {
    java -Xmx5000m \
      -jar ${PICARD} \
      CollectMultipleMetrics \
      INPUT=${in_bam} \
      REFERENCE_SEQUENCE=${ref_fa} \
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
    # If the queue system fails to move the stderr/stdout there must be something to cache:
    touch stderr
    touch stdout

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
  File in_bai
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_idx
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File ref_dict
  File ref_fa
  File ref_idx
  Int cpu=1
  File GATK
  String? U_option  # In case a -U option needs to be provided

  command {
    # Write the intervals to a file for better command debugging
    rand=`shuf -i 1-10000000 -n 1`
    mv ${write_lines(sequence_group_interval)} $rand.intervals

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx4000m \
      -jar ${GATK} \
      -T BaseRecalibrator \
      -R ${ref_fa} \
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
  File ref_fa
  File ref_idx
  Int cpu=1
  File GATK4

  command {
    # Write the intervals to a file for better command debugging
    rand=`shuf -i 1-10000000 -n 1`
    mv ${write_lines(sequence_group_interval)} $rand.intervals

    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
      -jar ${GATK4} ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fa} \
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
  Array[File] in_bqsr_reports
  String out_report_filename
  Int cpu=1
  File GATK

  command {
    java -Xmx3000m \
      -cp ${GATK} org.broadinstitute.gatk.tools.GatherBqsrReports \
      I=${sep=' I=' in_bqsr_reports} \
      O=${out_report_filename}
    }
  runtime {
    cpu: cpu
  }
  output {
    File out_bqsr_report = "${out_report_filename}"
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
  File ref_fa
  File ref_idx
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
      REFERENCE_SEQUENCE=${ref_fa} \
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
  File ref_fa
  File ref_idx
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectWgsMetrics \
      INPUT=${in_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fa} \
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
  File ref_fa
  File ref_idx
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectRawWgsMetrics \
      INPUT=${in_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fa} \
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
  File contamination_sites_vcf_idx
  String out_prefix
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
    --out ${out_prefix} \
    --bam ${in_bam} \
    1>/dev/null

    ${PYTHON3} <<CODE
    import csv
    import sys
    with open('${out_prefix}.selfSM') as selfSM:
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
    File selfSM = "${out_prefix}.selfSM"
    File depthSM = "${out_prefix}.depthSM"
    File log = "${out_prefix}.log"

    # We would like to do the following, however:
    # The object is read as a string
    # explicit string->float coercion via float(), as shown below, is supported by Cromwell
    # the interim value cannot be stored as a string and then assigned to a float. Variables intialized in output cannot be dereferenced in output.
    # Float contamination = float(read_object(${out_prefix} + ".selfSM").FREEMIX) / 0.75

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
  File ref_fa
  File ref_idx
  Float? contamination
  Int cpu=28
  File GATK

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar ${GATK} \
      -T HaplotypeCaller \
      -R ${ref_fa} \
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
    File out_gvcf = "${gvcf_basename}.vcf.gz"
    File out_gvcf_idx = "${gvcf_basename}.vcf.gz.tbi"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  Array[File] in_vcfs
  Array[File] in_vcfs_idxes
  String out_vcf_name
  Int cpu=1
  File PICARD

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xmx2g \
      -jar ${PICARD} \
      MergeVcfs \
      INPUT=${sep=' INPUT=' in_vcfs} \
      OUTPUT=${out_vcf_name}
  }
  output {
    File out_vcf = "${out_vcf_name}"
    File out_vcf_idx = "${out_vcf_name}.tbi"
  }
  runtime {
    cpu: cpu
  }
}

# Validate a GVCF with -gvcf specific validation
task ValidateGVCF {
  File in_vcf
  File in_vcf_idx
  File ref_fa
  File ref_idx
  File ref_dict
  File dbSNP_vcf
  File dbSNP_vcf_idx
  File wgs_calling_interval_list
  Int cpu=1
  File GATK

  command {
    java -Xmx8g \
      -jar ${GATK} \
      -T ValidateVariants \
      -V ${in_vcf} \
      -R ${ref_fa} \
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
  File in_vcf
  File in_vcf_idx
  String metrics_basename
  File dbSNP_vcf
  File dbSNP_vcf_idx
  File ref_dict
  Int cpu=1
  File PICARD

  command {
    java -Xmx2000m \
      -jar ${PICARD} \
      CollectVariantCallingMetrics \
      INPUT=${in_vcf} \
      OUTPUT=${metrics_basename} \
      DBSNP=${dbSNP_vcf} \
      SEQUENCE_DICTIONARY=${ref_dict} \
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
  File ref_fa
  File ref_idx
  String out_basename
  Int cpu=1
  File SAMTOOLS
  File seq_cache_populate

  command <<<
    set -e
    set -o pipefail

    ${SAMTOOLS} view -C -T ${ref_fa} ${in_bam} | \
    tee ${out_basename}.cram | \
    md5sum | awk '{print $1}' > ${out_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    ${seq_cache_populate} -root ./ref/cache ${ref_fa}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    ${SAMTOOLS} index ${out_basename}.cram
    mv ${out_basename}.cram.crai ${out_basename}.crai
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File out_cram = "${out_basename}.cram"
    File out_cram_idx = "${out_basename}.crai"
    File out_cram_md5 = "${out_basename}.cram.md5"
  }
}

# Convert a CRAM file to BAM format
task CramToBam {
  File ref_fa
  File ref_idx
  File ref_dict
  File cram_file
  String out_basename
  Int cpu=1
  File SAMTOOLS

command <<<
  set -e
  set -o pipefail

  ${SAMTOOLS} view -h -T ${ref_fa} ${cram_file} |
  ${SAMTOOLS} view -b -o ${out_basename}.bam -
  ${SAMTOOLS} index -b ${out_basename}.bam
  mv ${out_basename}.bam.bai ${out_basename}.bai
  >>>
  runtime {
    cpu: cpu
  }
  output {
    File out_bam = "${out_basename}.bam"
    File out_bai = "${out_basename}.bai"
  }
}
#############################
### TASK DEFINITIONS ENDS ###
#############################





############################
### WORKFLOW DEFINITIONS ###
############################
workflow WES_normal_tumor_somatic_SNV_wf {

  File contamination_sites_vcf
  File contamination_sites_vcf_idx

  String sample_name
  String base_file_name_normal
  String base_file_name_tumor
  String final_gvcf_ext
  File rawdata_normal_fastqR1
  File rawdata_normal_fastqR2
  File rawdata_tumor_fastqR1
  File rawdata_tumor_fastqR2
  String unmapped_bam_suffix

  Array[File] scattered_calling_intervals
  File wgs_calling_interval_list

  File ref_fa
  File ref_idx
  File ref_dict
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  File dbSNP_vcf
  File dbSNP_vcf_idx
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File gnomad_exome_vcf
  File gnomad_exome_vcf_idx
  File transcript_intervals

  String recalibrated_bam_basename_normal = base_file_name_normal + ".aligned.duplicates_marked.recalibrated"
  String recalibrated_bam_basename_tumor = base_file_name_tumor + ".aligned.duplicates_marked.recalibrated"
  String gvcf_name_normal = base_file_name_normal + final_gvcf_ext
  String gvcf_name_tumor = base_file_name_tumor + final_gvcf_ext
#  String final_gvcf_name_normal = base_file_name_normal + "_recall" + final_gvcf_ext
#  String final_gvcf_name_tumor = base_file_name_tumor + "_recall" + final_gvcf_ext

  # Tools
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

  # Variant recalibration:
  Array[String] SNP_annotations
  Array[String] INDEL_annotations
  Array[Float] SNP_tranches
  Array[Float] INDEL_tranches
  Array[String] SNP_resources
  Array[String] INDEL_resources
  Array[File] resource_files
  Array[File] resource_indices
  Float SNP_filter_level
  Float INDEL_filter_level

  # Custom hacks:
  String bwa_commandline = bwa + " mem -K 100000000 -p -v 3 -t 28 -Y $bash_ref_fa"
  String sub_strip_path = "/.*/"
  String sub_strip_unmapped = unmapped_bam_suffix + "$"


#######################
### WES NORMAL PART ###
#######################
  # Decompress and split the files into chunks, return an array of .fastq files
  call UnzipAndSplit as UnzipAndSplit_normal {
      input:
        PIGZ=pigz,
        in_fastqR1 = rawdata_normal_fastqR1,
        in_fastqR2 = rawdata_normal_fastqR2,
        sample_name = sample_name + '_normal'
  }


  Array[Pair[File, File]] fastqR1R2_chunks_normal = zip(UnzipAndSplit_normal.R1_splits, UnzipAndSplit_normal.R2_splits)
  scatter (fastq_chunk_normal in fastqR1R2_chunks_normal) {

    call TrimReads as TrimReads_normal {
        input:
          TRIMMOMATIC=trimmomatic,
          in_fastqR1 = fastq_chunk_normal.left,
          in_fastqR2 = fastq_chunk_normal.right,
          basenameR1 = basename(fastq_chunk_normal.left, ".fastq"),
          basenameR2 = basename(fastq_chunk_normal.right, ".fastq")
    }

    call FastqToBam as FastqToBam_normal {
        input:
          PICARD=picard,
          sample_name = sample_name + '_normal',
          in_fastqR1 = TrimReads_normal.out_R1,
          in_fastqR2 = TrimReads_normal.out_R2
    }

    # QC the unmapped BAM
    call CollectQualityYieldMetrics as CollectQualityYieldMetrics_normal {
      input:
        PICARD=picard,
        in_bam = FastqToBam_normal.out_bam,
        metrics_filename = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmapped.quality_yield_metrics"
    }

    # Map reads to reference
    call SamToFastqAndBwaMem as SamToFastqAndBwaMem_normal {
      input:
        SAMTOOLS=samtools,
        PICARD=picard,
        in_bam = FastqToBam_normal.out_bam,
        bwa_commandline = bwa_commandline,
        out_bam_basename = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmerged",
        ref_fa = ref_fa,
        ref_idx = ref_idx,
        ref_dict = ref_dict,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa
     }

    # Merge original uBAM and BWA-aligned BAM
    call MergeBamAlignment as MergeBamAlignment_normal {
      input:
        PICARD=picard,
        unmapped_bam = FastqToBam_normal.out_bam,
        aligned_bam = SamToFastqAndBwaMem_normal.out_bam,
        out_bam_basename = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".aligned.unsorted",
        ref_fa = ref_fa,
        ref_idx = ref_idx,
        ref_dict = ref_dict
    }

    # QC the aligned but unsorted readgroup BAM
    # No reference needed as the input here is unsorted; providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics_normal {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment_normal.out_bam,
        out_bam_prefix = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".readgroup"
    }

    # Sort and fix tags in the merged BAM
    call SortAndFixTags as SortAndFixReadGroupBam_normal {
      input:
        PICARD=picard,
        in_bam = MergeBamAlignment_normal.out_bam,
        out_bam_basename = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".sorted",
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }

    # Validate the aligned and sorted readgroup BAM
    # This is called to help in finding problems early.
    # If considered too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile_normal {
      input:
        PICARD=picard,
        ref_fa = ref_fa,
        ref_idx = ref_idx,
        ref_dict = ref_dict,
        in_bam = SortAndFixReadGroupBam_normal.out_bam,
        in_bai = SortAndFixReadGroupBam_normal.out_bai,
        report_filename = sub(sub(FastqToBam_normal.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".validation_report"
    }

  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates as MarkDuplicates_normal {
    input:
      PICARD=picard,
      in_bams = MergeBamAlignment_normal.out_bam,
      out_bam_basename = base_file_name_normal + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name_normal + ".duplicate_metrics"
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam_normal {
    input:
      PICARD=picard,
      in_bam = MarkDuplicates_normal.out_bam,
      out_bam_basename = base_file_name_normal + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV as CreateSequenceGroupingTSV_normal {
    input:
      ref_dict = ref_dict
  }

  # Estimate level of cross-sample contamination
  call CheckContamination as CheckContamination_normal {
    input:
      verifyBamID=verifyBamID,
      PYTHON3 = python3,
      in_bam = SortAndFixSampleBam_normal.out_bam,
      in_bai = SortAndFixSampleBam_normal.out_bai,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_idx = contamination_sites_vcf_idx,
      out_prefix = base_file_name_normal + ".preBqsr"
  }

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV_normal.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator as BaseRecalibrator_normal {
      input:
        GATK=gatk,
        in_bam = SortAndFixSampleBam_normal.out_bam,
        in_bai = SortAndFixSampleBam_normal.out_bai,
        recalibration_report_filename = base_file_name_normal + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_idx = dbSNP_vcf_idx,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports as GatherBqsrReports_normal {
    input:
      GATK=gatk,
      in_bqsr_reports = BaseRecalibrator_normal.recalibration_report,
      out_report_filename = base_file_name_normal + ".recal_data.csv"
  }

  scatter (subgroup in CreateSequenceGroupingTSV_normal.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR as ApplyBQSR_normal {
      input:
        GATK4=gatk4,
        in_bam = SortAndFixSampleBam_normal.out_bam,
        in_bai = SortAndFixSampleBam_normal.out_bai,
        out_bam_basename = recalibrated_bam_basename_normal,
        recalibration_report = GatherBqsrReports_normal.out_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles as GatherBamFiles_normal {
    input:
      PICARD=picard,
      in_bams = ApplyBQSR_normal.recalibrated_bam,
      out_bam_basename = base_file_name_normal
  }

  # QC the final BAM (consolidated after scattered BQSR)
  call CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      out_bam_prefix = base_file_name_normal + ".readgroup",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # Validate the final BAM
  call ValidateSamFile as ValidateAggregatedSamFile_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      report_filename = base_file_name_normal + ".validation_report",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics as CollectAggregationMetrics_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      out_bam_prefix = base_file_name_normal,
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # QC the sample WGS metrics (stringent thresholds)
  call CollectWgsMetrics as CollectWgsMetrics_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      metrics_filename = base_file_name_normal + ".wgs_metrics",
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # QC the sample raw WGS metrics (common thresholds)
  call CollectRawWgsMetrics as CollectRawWgsMetrics_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      metrics_filename = base_file_name_normal + ".raw_wgs_metrics",
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # Generate a checksum per readgroup in the final BAM
  call CalculateReadGroupChecksum as CalculateReadGroupChecksum_normal {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_normal.out_bam,
      in_bai = GatherBamFiles_normal.out_bai,
      read_group_md5_filename = recalibrated_bam_basename_normal + ".bam.read_group_md5"
  }

  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram as ConvertToCram_normal {
    input:
      seq_cache_populate=seq_cache_populate,
      SAMTOOLS=samtools,
      in_bam = GatherBamFiles_normal.out_bam,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      out_basename = base_file_name_normal
  }

  # Convert the CRAM back to BAM to check that the conversions do not introduce errors
  call CramToBam as CramToBam_normal {
    input:
      SAMTOOLS=samtools,
      ref_fa = ref_fa,
      ref_dict = ref_dict,
      ref_idx = ref_idx,
      cram_file = ConvertToCram_normal.out_cram,
      out_basename = base_file_name_normal + ".roundtrip"
  }

  # Validate the roundtripped BAM
  call ValidateSamFile as ValidateBamFromCram_normal {
    input:
      PICARD=picard,
      in_bam = CramToBam_normal.out_bam,
      in_bai = CramToBam_normal.out_bai,
      report_filename = base_file_name_normal + ".bam.roundtrip.validation_report",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      max_output = 1000000000,
      ignore = ["null"]
  }

  # Call variants in parallel over WGS calling intervals
  scatter (subInterval in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller as HaplotypeCaller_normal {
      input:
        GATK=gatk,
        contamination = CheckContamination_normal.contamination,
        in_bam = GatherBamFiles_normal.out_bam,
        in_bai = GatherBamFiles_normal.out_bai,
        interval_list = subInterval,
        gvcf_basename = base_file_name_normal,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
     }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs as MergeVCFs_normal {
    input:
      PICARD=picard,
      in_vcfs = HaplotypeCaller_normal.out_gvcf,
      in_vcfs_idxes = HaplotypeCaller_normal.out_gvcf_idx,
      out_vcf_name = gvcf_name_normal
  }

  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF as ValidateGVCF_normal {
    input:
      GATK=gatk,
      in_vcf = MergeVCFs_normal.out_vcf,
      in_vcf_idx = MergeVCFs_normal.out_vcf_idx,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list
  }

  # QC the GVCF
  call CollectGvcfCallingMetrics as CollectGvcfCallingMetrics_normal {
    input:
      PICARD=picard,
      in_vcf = MergeVCFs_normal.out_vcf,
      in_vcf_idx = MergeVCFs_normal.out_vcf_idx,
      metrics_basename = base_file_name_normal,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx,
      ref_dict = ref_dict
  }



#######################
# Experiemental VQSR
############################

  # Build SNP model:
  call BuildVQSRModel as BuildVQSRModelForSNPs_normal {
    input:
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      in_vcf = MergeVCFs_normal.out_vcf,
      in_vcf_idx = MergeVCFs_normal.out_vcf_idx,
      wgs_calling_interval_list = wgs_calling_interval_list,
      out_basename = base_file_name_normal,
      annotations = SNP_annotations,
      mode = "SNP",
      tranches = SNP_tranches,
      resources = SNP_resources,
      resource_files = resource_files,
      resource_indices = resource_indices
  }

  # Build INDEL model:
  call BuildVQSRModel as BuildVQSRModelForINDELs_normal {
    input:
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      in_vcf = MergeVCFs_normal.out_vcf,
      in_vcf_idx = MergeVCFs_normal.out_vcf_idx,
      wgs_calling_interval_list = wgs_calling_interval_list,
      out_basename = base_file_name_normal,
      annotations = INDEL_annotations,
      mode = "INDEL",
      tranches = INDEL_tranches,
      resources = INDEL_resources,
      resource_files = resource_files,
      resource_indices = resource_indices
  }

  # Apply INDEL filter (first because INDEL model is usually done sooner):
  call ApplyRecalibrationFilter as ApplyRecalibrationFilterForINDELs_normal {
    input:
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      in_vcf = MergeVCFs_normal.out_vcf,
      in_vcf_idx = MergeVCFs_normal.out_vcf_idx,
      wgs_calling_interval_list = wgs_calling_interval_list,
      out_basename = base_file_name_normal + ".recal.INDEL",
      mode = "INDEL",
      recal_file = BuildVQSRModelForINDELs.recal_file,
      recal_file_idx = BuildVQSRModelForINDELs.recal_file_idx,
      tranches_file = BuildVQSRModelForINDELs.tranches_file,
      filter_level = INDEL_filter_level
  }

  # Apply SNP filter:
  call ApplyRecalibrationFilter as ApplyRecalibrationFilterForSNPs_normal {
    input:
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      in_vcf = ApplyRecalibrationFilterForINDELs.out_vcf,
      in_vcf_idx = ApplyRecalibrationFilterForINDELs.out_vcf_idx,
      wgs_calling_interval_list = wgs_calling_interval_list,
      out_basename = base_file_name_normal + ".recal.final",
      mode = "SNP",
      recal_file = BuildVQSRModelForSNPs.recal_file,
      recal_file_idx = BuildVQSRModelForSNPs.recal_file_idx,
      tranches_file = BuildVQSRModelForSNPs.tranches_file,
      filter_level = SNP_filter_level
  }


##############################################################################



############################
### WES NORMAL PART ENDS ###
############################







######################
### WES TUMOR PART ###
######################
  # Decompress and split the files into chunks, return an array of .fastq files
  call UnzipAndSplit as UnzipAndSplit_tumor {
      input:
        PIGZ=pigz,
        in_fastqR1 = rawdata_tumor_fastqR1,
        in_fastqR2 = rawdata_tumor_fastqR2,
        sample_name = sample_name + '_tumor'
  }

  Array[Pair[File, File]] fastqR1R2_chunks_tumor = zip(UnzipAndSplit_tumor.R1_splits, UnzipAndSplit_tumor.R2_splits)
  scatter (fastq_chunk_tumor in fastqR1R2_chunks_tumor) {
    call TrimReads as TrimReads_tumor {
        input:
          TRIMMOMATIC=trimmomatic,
          in_fastqR1 = fastq_chunk_tumor.left,
          in_fastqR2 = fastq_chunk_tumor.right,
          basenameR1 = basename(fastq_chunk_tumor.left, ".fastq"),
          basenameR2 = basename(fastq_chunk_tumor.right, ".fastq")
    }

    # Convert to unmapped BAM file
    call FastqToBam as FastqToBam_tumor {
        input:
          PICARD=picard,
          sample_name = sample_name + '_tumor',
          in_fastqR1 = TrimReads_tumor.out_R1,
          in_fastqR2 = TrimReads_tumor.out_R2
    }

    # QC the unmapped BAM
    call CollectQualityYieldMetrics as CollectQualityYieldMetrics_tumor {
      input:
        PICARD=picard,
        in_bam = FastqToBam_tumor.out_bam,
        metrics_filename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmapped.quality_yield_metrics"
    }

    # Map reads to reference
    call SamToFastqAndBwaMem as SamToFastqAndBwaMem_tumor {
      input:
        SAMTOOLS=samtools,
        PICARD=picard,
        in_bam = FastqToBam_tumor.out_bam,
        bwa_commandline = bwa_commandline,
        out_bam_basename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmerged",
        ref_fa = ref_fa,
        ref_idx = ref_idx,
        ref_dict = ref_dict,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa
     }

    # Merge original uBAM and BWA-aligned BAM
    call MergeBamAlignment as MergeBamAlignment_tumor {
      input:
        PICARD=picard,
        unmapped_bam = FastqToBam_tumor.out_bam,
        aligned_bam = SamToFastqAndBwaMem_tumor.out_bam,
        out_bam_basename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".aligned.unsorted",
        ref_fa = ref_fa,
        ref_idx = ref_idx,
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
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }

    # Validate the aligned and sorted readgroup BAM
    # This is called to help in finding problems early.
    # If considered too time consuming and not helpful, can be removed.
    call ValidateSamFile as ValidateReadGroupSamFile_tumor {
      input:
        PICARD=picard,
        ref_fa = ref_fa,
        ref_idx = ref_idx,
        ref_dict = ref_dict,
        in_bam = SortAndFixReadGroupBam_tumor.out_bam,
        in_bai = SortAndFixReadGroupBam_tumor.out_bai,
        report_filename = sub(sub(FastqToBam_tumor.out_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".validation_report"
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

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags as SortAndFixSampleBam_tumor {
    input:
      PICARD=picard,
      in_bam = MarkDuplicates_tumor.out_bam,
      out_bam_basename = base_file_name_tumor + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  call CreateSequenceGroupingTSV as CreateSequenceGroupingTSV_tumor {
    input:
      ref_dict = ref_dict
  }

  # Estimate level of cross-sample contamination
  call CheckContamination as CheckContamination_tumor {
    input:
      verifyBamID=verifyBamID,
      PYTHON3 = python3,
      in_bam = SortAndFixSampleBam_tumor.out_bam,
      in_bai = SortAndFixSampleBam_tumor.out_bai,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_idx = contamination_sites_vcf_idx,
      out_prefix = base_file_name_tumor + ".preBqsr"
  }

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV_tumor.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator as BaseRecalibrator_tumor {
      input:
        GATK=gatk,
        in_bam = SortAndFixSampleBam_tumor.out_bam,
        in_bai = SortAndFixSampleBam_tumor.out_bai,
        recalibration_report_filename = base_file_name_tumor + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_idx = dbSNP_vcf_idx,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  call GatherBqsrReports as GatherBqsrReports_tumor {
    input:
      GATK=gatk,
      in_bqsr_reports = BaseRecalibrator_tumor.recalibration_report,
      out_report_filename = base_file_name_tumor + ".recal_data.csv"
  }

  scatter (subgroup in CreateSequenceGroupingTSV_tumor.sequence_grouping_with_unmapped) {

    # Apply the recalibration model by interval
    call ApplyBQSR as ApplyBQSR_tumor {
      input:
        GATK4=gatk4,
        in_bam = SortAndFixSampleBam_tumor.out_bam,
        in_bai = SortAndFixSampleBam_tumor.out_bai,
        out_bam_basename = recalibrated_bam_basename_tumor,
        recalibration_report = GatherBqsrReports_tumor.out_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles as GatherBamFiles_tumor {
    input:
      PICARD=picard,
      in_bams = ApplyBQSR_tumor.recalibrated_bam,
      out_bam_basename = base_file_name_tumor
  }

  # QC the final BAM (consolidated after scattered BQSR)
  call CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      out_bam_prefix = base_file_name_tumor + ".readgroup",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # Validate the final BAM
  call ValidateSamFile as ValidateAggregatedSamFile_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      report_filename = base_file_name_tumor + ".validation_report",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      ignore = ["null"]
  }

  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics as CollectAggregationMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      out_bam_prefix = base_file_name_tumor,
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # QC the sample WGS metrics (stringent thresholds)
  call CollectWgsMetrics as CollectWgsMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      metrics_filename = base_file_name_tumor + ".wgs_metrics",
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # QC the sample raw WGS metrics (common thresholds)
  call CollectRawWgsMetrics as CollectRawWgsMetrics_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      metrics_filename = base_file_name_tumor + ".raw_wgs_metrics",
      ref_fa = ref_fa,
      ref_idx = ref_idx
  }

  # Generate a checksum per readgroup in the final BAM
  call CalculateReadGroupChecksum as CalculateReadGroupChecksum_tumor {
    input:
      PICARD=picard,
      in_bam = GatherBamFiles_tumor.out_bam,
      in_bai = GatherBamFiles_tumor.out_bai,
      read_group_md5_filename = recalibrated_bam_basename_tumor + ".bam.read_group_md5"
  }

  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram as ConvertToCram_tumor {
    input:
      seq_cache_populate=seq_cache_populate,
      SAMTOOLS=samtools,
      in_bam = GatherBamFiles_tumor.out_bam,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      out_basename = base_file_name_tumor
  }

  # Convert the CRAM back to BAM to check that the conversions do not introduce errors
  call CramToBam as CramToBam_tumor {
    input:
      SAMTOOLS=samtools,
      ref_fa = ref_fa,
      ref_dict = ref_dict,
      ref_idx = ref_idx,
      cram_file = ConvertToCram_tumor.out_cram,
      out_basename = base_file_name_tumor + ".roundtrip"
  }

  # Validate the roundtripped BAM
  call ValidateSamFile as ValidateBamFromCram_tumor {
    input:
      PICARD=picard,
      in_bam = CramToBam_tumor.out_bam,
      in_bai = CramToBam_tumor.out_bai,
      report_filename = base_file_name_tumor + ".bam.roundtrip.validation_report",
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      max_output = 1000000000,
      ignore = ["null"]
  }

  # Call variants in parallel over WGS calling intervals
  scatter (subInterval in scattered_calling_intervals) {

    # Generate GVCF by interval
    call HaplotypeCaller as HaplotypeCaller_tumor {
      input:
        GATK=gatk,
        contamination = CheckContamination_tumor.contamination,
        in_bam = GatherBamFiles_tumor.out_bam,
        in_bai = GatherBamFiles_tumor.out_bai,
        interval_list = subInterval,
        gvcf_basename = base_file_name_tumor,
        ref_dict = ref_dict,
        ref_fa = ref_fa,
        ref_idx = ref_idx
     }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call MergeVCFs as MergeVCFs_tumor {
    input:
      PICARD=picard,
      in_vcfs = HaplotypeCaller_tumor.out_gvcf,
      in_vcfs_idxes = HaplotypeCaller_tumor.out_gvcf_idx,
      out_vcf_name = gvcf_name_tumor
  }

  # Validate the GVCF output of HaplotypeCaller
  call ValidateGVCF as ValidateGVCF_tumor {
    input:
      GATK=gatk,
      in_vcf = MergeVCFs_tumor.out_vcf,
      in_vcf_idx = MergeVCFs_tumor.out_vcf_idx,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list
  }

  # QC the GVCF
  call CollectGvcfCallingMetrics as CollectGvcfCallingMetrics_tumor {
    input:
      PICARD=picard,
      in_vcf = MergeVCFs_tumor.out_vcf,
      in_vcf_idx = MergeVCFs_tumor.out_vcf_idx,
      metrics_basename = base_file_name_tumor,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx,
      ref_dict = ref_dict,
  }

###########################
### WES TUMOR PART ENDS ###
###########################





########################
### SOMATIC VARIANTS ###
########################

  # Somatic variant calling with MuTect2
  call MuTect2 {
    input:
      GATK4_LAUNCH=gatk4_launch,
      contamination = CheckContamination_tumor.contamination,
      in_bam_tumor = GatherBamFiles_tumor.out_bam,
      in_bai_tumor = GatherBamFiles_tumor.out_bai,
      in_bam_normal = GatherBamFiles_normal.out_bam,
      in_bai_normal = GatherBamFiles_normal.out_bai,
      sample_name = sample_name,
      sample_name_tumor = sample_name+'_tumor',
      sample_name_normal = sample_name+'_normal',
      ref_dict = ref_dict,
      ref_fa = ref_fa,
      ref_idx = ref_idx,
      transcript_intervals = transcript_intervals,
      gnomad_exome_vcf = gnomad_exome_vcf,
      gnomad_exome_vcf_idx = gnomad_exome_vcf_idx,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx
   }

  # Filter somatic variants
  call FilterMutectCalls {
    input:
      GATK4_LAUNCH=gatk4_launch,
      contamination = CheckContamination_tumor.contamination,
      sample_name = sample_name,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_idx = dbSNP_vcf_idx,
      in_vcf = MuTect2.out_vcf,
      in_vcf_idx = MuTect2.out_vcf_idx
   }

  # Filter somatic variants second round
  call FilterByOrientationBias {
    input:
      GATK4_LAUNCH=gatk4_launch,
      sample_name = sample_name,
      in_vcf = FilterMutectCalls.out_vcf,
      in_vcf_idx = FilterMutectCalls.out_vcf_idx,
      pre_adapter_detail_metrics_tumor = CollectAggregationMetrics_tumor.pre_adapter_detail_metrics
   }

#############################
### SOMATIC VARIANTS ENDS ###
#############################



  # Outputs that will be retained when execution is complete
  output {
    CollectQualityYieldMetrics_normal.*
    CollectQualityYieldMetrics_tumor.*
    ValidateReadGroupSamFile_normal.*
    ValidateReadGroupSamFile_tumor.*
    CollectReadgroupBamQualityMetrics_normal.*
    CollectReadgroupBamQualityMetrics_tumor.*
    CollectUnsortedReadgroupBamQualityMetrics_normal.*
    CollectUnsortedReadgroupBamQualityMetrics_tumor.*
    ValidateBamFromCram_normal.*
    ValidateBamFromCram_tumor.*
    CalculateReadGroupChecksum_normal.*
    CalculateReadGroupChecksum_tumor.*
    ValidateAggregatedSamFile_normal.*
    ValidateAggregatedSamFile_tumor.*
    CollectAggregationMetrics_normal.*
    CollectAggregationMetrics_tumor.*
    CollectWgsMetrics_normal.*
    CollectWgsMetrics_tumor.*
    CollectRawWgsMetrics_normal.*
    CollectRawWgsMetrics_tumor.*
    CheckContamination_normal.*
    CheckContamination_tumor.*
    CollectGvcfCallingMetrics_normal.*
    CollectGvcfCallingMetrics_tumor.*
    MarkDuplicates_normal.duplicate_metrics
    MarkDuplicates_tumor.duplicate_metrics
    GatherBqsrReports_normal.*
    GatherBqsrReports_tumor.*
    ConvertToCram_normal.*
    ConvertToCram_tumor.*
    MergeVCFs_normal.*
    MergeVCFs_tumor.*
    FilterByOrientationBias.*
    BuildVQSRModelForSNPs_normal.*
    BuildVQSRModelForINDELs_normal.*
    ApplyRecalibrationFilterForINDELs_normal.*
  }
}
#################################
### WORKFLOW DEFINITIONS ENDS ###
#################################
