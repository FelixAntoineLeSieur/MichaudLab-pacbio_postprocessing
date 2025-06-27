//This process is based from the Nf-core bcftools_isec module:
//https://nf-co.re/modules/bcftools_isec

//It was used as a local module to make some modifications, mostly inclusion of the final report 
process VALIDATION_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("${prefix}", type: "dir"), emit: results
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    bcftools isec  \\
        $args \\
        -p $prefix \\
        ${vcfs}

    ##Report
    siteFile=$prefix/sites.txt
    echo "Report for samplename" >"${prefix}/${prefix}_summary.txt"
    uniqList=`cut \$siteFile -f 1-2 | uniq -c`
    numVar=`echo "\$uniqList" | wc -l`

    echo "Num of Non-Ref variants among 44 positions:\$numVar" >>"${prefix}/${prefix}_summary.txt"
    mismatchList=`echo \$uniqList | grep "2 "`
    mismatchListNoC=`echo "\$mismatchList"| cut -d' ' -f8-`
    misfitCount=`printf "\$mismatchList"| wc -l`
    echo "Num of variants of same position not matching :\$misfitCount" >>"${prefix}/${prefix}_summary.txt"
    if [ \$misfitCount != 0 ]; then
        printf "Misfit list: \n\$mismatchListNoC\n" >"${prefix}/${prefix}_report.txt"
    else
        printf "No mismatched genotypes to report\n" >"${prefix}/${prefix}_report.txt"
    fi
    nonMisfit=`grep -v "\$mismatchListNoC" \$siteFile`

    longreadOnly=`echo "\$nonMisfit" | awk '/\t10\$/{print \$0}'`
    longreadCoords=`echo "\$longreadOnly" | cut -f 1-4`
    longreadOnlyCounts=`printf "\$longreadCoords" | wc -l`
    echo "Num of variants only detected in long reads :\$longreadOnlyCounts" >>"${prefix}/${prefix}_summary.txt"
    printf "LongReads only list: \n\$longreadOnly\n" >>"${prefix}/${prefix}_report.txt"
    
    shortreadOnly=`echo "\$nonMisfit" | awk '/\t01\$/{print \$0}'`
    shortreadOnlyCounts=`echo "\$shortreadOnly" | wc -l`
    shortreadCoords=`echo "\$shortreadOnly" | cut -f 1-4`
    shortreadOnlyCounts=`printf "\$shortreadCoords" | wc -l`
    echo "Num of variants only detected in short reads :\$shortreadOnlyCounts" >>"${prefix}/${prefix}_summary.txt"
    printf "ShortReads only list: \n\$shortreadOnly\n" >>"${prefix}/${prefix}_report.txt"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS

    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/README.txt
    touch ${prefix}/sites.txt
    echo "" | gzip > ${prefix}/0000.vcf.gz
    touch ${prefix}/0000.vcf.gz.tbi
    echo "" | gzip > ${prefix}/0001.vcf.gz
    touch ${prefix}/0001.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
