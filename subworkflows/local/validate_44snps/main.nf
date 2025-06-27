// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

// include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_LR } from '../../../modules/nf-core/gatk4/selectvariants/main' 
include { GATK4_SELECTVARIANTS as GATK4_SELECTVARIANTS_SR } from '../../../modules/nf-core/gatk4/selectvariants/main' 
include { GATK4_GENOTYPEGVCFS }                             from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { VALIDATION_REPORT }                               from '../../../modules/local/validationreport/main'

workflow VALIDATE_44SNPS {

    take:
    ch_validate_lr // channel: [ val(meta), [ lr_VCF ] [lr_VCF_index] [IntervalBED] ]
    ch_validate_sr
    genomeFasta
    genomeFai
    genomeDict

    main:

    ch_versions = Channel.empty()
    ch_select_out_lr = GATK4_SELECTVARIANTS_LR(ch_validate_lr).vcf.join(GATK4_SELECTVARIANTS_LR.out.tbi)
    ch_select_out_sr = GATK4_SELECTVARIANTS_SR(ch_validate_sr).vcf.join(GATK4_SELECTVARIANTS_SR.out.tbi)

    ch_input_genotypegvcf_sr = ch_select_out_sr.map{meta,vcf,tbi -> [meta,vcf,tbi, [], []]}
    ch_output_genotypegvcf_sr = GATK4_GENOTYPEGVCFS(
        ch_input_genotypegvcf_sr,
        [[:], genomeFasta],
        [[:], genomeFai],
        [[:], genomeDict],
        [[:], []],
        [[:], []]
        ).vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
    
    ch_input_isec = ch_select_out_lr.join(ch_output_genotypegvcf_sr).map {
        meta,lrVcf,lrIndex,srVcf,srIndex ->
            [meta,[lrVcf,srVcf],[lrIndex,srIndex]]
    }
    VALIDATION_REPORT(ch_input_isec)

    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS_LR.out.versions.first())
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())
    ch_versions = ch_versions.mix(VALIDATION_REPORT.out.versions.first())
    // SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    // ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
