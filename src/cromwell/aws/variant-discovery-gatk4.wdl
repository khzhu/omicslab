version 1.0

## Copyright Arcus Omics Team, 2021
## 
## This WDL pipeline implements germline variant discovery according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
## - Germline variants in a VCF format and annovar annotated variants in a TSV format.
##
## Software version requirements 
## - GATK 4 or later
## - BWA 0.7.15-r1140
## - Picard 2.16.0-SNAPSHOT
## - Samtools 1.3.1 (using htslib 1.3.1)
## - Annovar 2018.04.18
##
## Cromwell version support 
## - Successfully tested on v59
##
## Runtime parameters are optimized for Amazon's Genomic Cli Platform implementation.
##

import "imports/bwa-map2genome.wdl" as bwa_mapper
import "imports/bam2cram.wdl" as utils
import "imports/haplotypecaller-gvcf.wdl" as haplo_caller
import "imports/annovar-consensus.wdl" as anno_consensus

# WORKFLOW DEFINITION 
workflow VariantDiscovery_GATK4 {
  input {
    String sample_name
    String ref_name
    String unmapped_bam_suffix
    File unmapped_bam
    File scattered_calling_intervals_list
    #Human Genome
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_sa 
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_amb
    #dbSNP
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    #Annotation
    File annovar_tar_ball
    String annovar_protocols
    String annovar_operation
  }

  # Map reads to reference
  call bwa_mapper.BWAMap2GenomeWF {
    input:
      sample_name = sample_name,
      ref_name = ref_name,
      unmapped_bam_suffix = unmapped_bam_suffix,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_sa = ref_sa,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_pac = ref_pac,
      ref_amb = ref_amb,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      unmapped_bam = unmapped_bam
  }

  ## convert a bwamem mapped bam file to a cram
  call utils.BamToCramWF {
    input:
      input_bam = BWAMap2GenomeWF.analysis_ready_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      base_file_name = sample_name + "." + ref_name
  }

  ## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
  ## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
  call haplo_caller.Haplotype2GVCF {
    input:
      input_bam = BWAMap2GenomeWF.analysis_ready_bam,
      input_bam_index = BWAMap2GenomeWF.analysis_ready_bam_index,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      scattered_calling_intervals_list = scattered_calling_intervals_list
  }

  ## Consensus variants for human DNA sequencing using annovar software
  call anno_consensus.AnnovarConsensusWF {
    input:
        ref_name = ref_name,
        gvcf_file = Haplotype2GVCF.output_vcf,
        annovar_tar_ball = annovar_tar_ball,
        annovar_protocols = annovar_protocols,
        annovar_operation = annovar_operation
  }

  ## Outputs that will be retained when execution is complete  
  output {
    File output_GATK_annotated_vcf = AnnovarConsensusWF.GATK_annotated_vcf
    File output_GATK_annotated_table = AnnovarConsensusWF.GATK_annotated
    File output_gvcf_file = Haplotype2GVCF.output_vcf
    File output_gvcf_file_index = Haplotype2GVCF.output_vcf_index
    File output_cram = BamToCramWF.output_cram
    File output_cram_index = BamToCramWF.output_cram_index
    File output_cram_md5 = BamToCramWF.output_cram_md5
  }
}
