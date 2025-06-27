/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
import groovy.json.JsonSlurper

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PEDDY                  } from '../modules/nf-core/peddy/main'
include { VALIDATE_44SNPS        } from '../subworkflows/local/validate_44snps'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pacbio_postprocessing_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PACBIO_POSTPROCESSING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    reports          = Channel.empty()



    def genomeFasta = file(params.genomeFasta)
    def genomeFai = file(genomeFasta + ".fai")
    def genomeDict = file(params.genomeFasta.substring(0,params.genomeFasta.indexOf(".")) + ".dict")

    //Subworkflow to validate our SNV vcf with the short reads GVCF
    if (params.validate44SNP && file(params.validationSNPBed, checkIfExists: true).exists()){
        ch_validate_input = ch_samplesheet.multiMap{
            meta,output,rawBam,gvcf,bam_haplo ->
            def jsonSlurper = new JsonSlurper()
            def outputJSON = jsonSlurper.parseText(file(output,checkIfExists: true).text)
            //Old format of outputs json
            if (outputJSON.humanwgs_sample_phased_small_variant_vcfs){
                lrVCF = outputJSON.humanwgs_sample_phased_small_variant_vcfs.data
                lrVCF_index = outputJSON.humanwgs_sample_phased_small_variant_vcfs.data_index
            } else if (outputJSON.humanwgs_singleton_phased_small_variant_vcf){
                //Newer format distinguishes between family and singleton
                if ((meta.outputPrefix == 'humanwgs_singleton') && file(outputJSON.humanwgs_singleton_phased_small_variant_vcf,checkIfExists: true).exists()){
                lrVCF = file(outputJSON.humanwgs_singleton_phased_small_variant_vcf,checkIfExists: true)
                lrVCF_index = file(outputJSON.humanwgs_singleton_phased_small_variant_vcf_index,checkIfExists: true)
                } else if ((meta.outputPrefix == 'humanwgs_family') && file(outputJSON.humanwgs_family_phased_small_variant_vcf,checkIfExists: true).exists()){
                file(outputJSON.humanwgs_family_phased_small_variant_vcf,checkIfExists: true)
                lrVCF = file(outputJSON.humanwgs_singleton_phased_small_variant_vcf,checkIfExists: true)
                lrVCF_index = file(outputJSON.humanwgs_singleton_phased_small_variant_vcf_index,checkIfExists: true)
                }
            else{
                error("The VCF file was not found in the output json for sample $meta.id")
            }
            }
            srGVCF = gvcf
            srGVCF_index = ((srGVCF != []) && (file(srGVCF + ".tbi").exists())) ? file(srGVCF + ".tbi") : []
            lr:[meta,[lrVCF],[lrVCF_index],[file(params.validationSNPBed,checkIfExists: true)]]
            sr:[meta,[srGVCF],[srGVCF_index],[file(params.validationSNPBed,checkIfExists: true)]]
        }

        ch_validate_output = VALIDATE_44SNPS(ch_validate_input,genomeFasta,genomeFai,genomeDict)
    }

    if (params.inputFamily) {
    // PEDDY module
    chPEDDYInput = ch_samplesheet.map{
        meta, output, bam -> 
        [meta,output.humanwgs_family_pedigree]}
        
    }

    //Gather all possible reports for MULTIQC

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pacbio_postprocessing_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }



    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    ch_multiqc_files = ch_multiqc_files.mix(ch_samplesheet.map{meta,output,rawBam,gvcf,bam_haplo -> 
        [output.parent.parent]}).collect()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
