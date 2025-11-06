#!/usr/bin/env nextflow
workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}


workflow.onError = {
    println "Oops .. something went wrong"
}

params.help=false
params.in = ""
params.out = ""
params.cpus =  4
params.fast = false
params.check_catalogue = false
params.catalogue_name = ""
params.catalogue = ""
params.minimum_memory = 30.GB
params.allowed_catalogues = "fc_1_3_gut,gg_13_6_caecal,clf_1_0_gut,hs_10_4_gut,hs_8_4_oral,hs_2_9_skin,mm_5_0_gut,oc_5_7_gut,rn_5_9_gut,ssc_9_3_gut"

def usage() {
    println("nf-meteor.nf --in <fastq_dir> --catalogue_name <catalogue_name> --out <output_dir> --cpus <nb_cpus> -w <temp_work_dir>")
    println("--in Directory containing paired fastq.gz files (default ${params.in}).")
    println("--out Output directory (default ${params.out}). ")
    println("--cpus Number of cpus to use (default ${params.cpus}).")
    println("--catalogue_name Name of the prebuilt catalogue to use (default: none). Allowed values are: ${params.allowed_catalogues}")
    println("--catalogue Path to a custom catalogue (overrides --catalogue_name if both are provided).")
    println("--fast Enable fast mode for meteor (no functional analysis) (default: ${params.fast}).")
    println("--check_catalogue Check md5sum of the catalogue is compatible with the input reads (default: ${params.check_catalogue}).")
}

// Convert string to list for validation
def allowed_catalogues = params.allowed_catalogues.split(',')

if(params.help){
    usage()
    exit(1)
}

// Validate catalogue_name parameter
if (params.catalogue_name && !allowed_catalogues.contains(params.catalogue_name)) {
    println "ERROR: Invalid catalogue name '${params.catalogue_name}'"
    println "Allowed catalogues are:"
    allowed_catalogues.each { println "  - ${it}" }
    exit 1
}


myDir = file(params.out)
myDir.mkdirs()

process meteor_download {
    tag { params.catalogue_name }
    conda "meteor=2.0.21"
    
    output:
    path("${params.catalogue_name}${params.fast ? '_taxo' : ''}"), emit: catalogue

    script:
    def fast_option = params.fast ? "--fast" : ""
    def check_catalogue = params.check_catalogue ? "-c" : ""
    """
    meteor download -i ${params.catalogue_name} -o ./ ${check_catalogue} ${fast_option}
    """
}

process meteor_fastq {
    tag { reads_id }
    conda "meteor=2.0.21"

    input:
    tuple val(reads_id), path(forward), path(reverse), val(count)

    output:
    tuple val(reads_id), path("fastq/*"), val(count)

    """
    meteor fastq -i ./ -p -o fastq
    """
}

process meteor_mapping {
    tag { reads_id }
    cpus params.cpus
    memory { 
        def baseMemoryGB = params.minimum_memory.toGiga()
        def totalMemoryGB = baseMemoryGB + (count as Double)
        
        if (params.fast) {
            // Fast mode: maximum 10G
            "${Math.min(totalMemoryGB, 10)}G"
        } else {
            // Normal mode: original calculation
            "${totalMemoryGB}G"
        }
    }
    conda "meteor=2.0.21"

    input:
    tuple val(reads_id), path(fastq), val(count)
    path(catalogue)

    output:
    tuple val(reads_id), path("mapping/*"), emit: mapping

    """
    meteor mapping -i ${fastq} -r ${catalogue} -t ${params.cpus} -o mapping --kf
    """
}

process meteor_profile {
    tag { reads_id }
    cpus params.cpus
    conda "meteor=2.0.21"

    input:
    tuple val(reads_id), path(mapping)
    path(catalogue)

    output:
    tuple val(reads_id), path("profile/*"), emit: profile

    """
    meteor profile -i ${mapping} -r ${catalogue} -o profile
    """
}

process meteor_merge {
    memory { 
        def samples = profile instanceof List ? profile.size() : 1
        // 200MB per sample with 1.5x safety margin
        def memoryGB = Math.ceil(samples * 0.2 * 1.5) 
        return "${memoryGB}G"
    }
    conda "meteor=2.0.21"
    publishDir "$myDir", mode: 'copy'

    input:
    path(profile)
    path(catalogue)

    output:
    path("merged")

    """
    meteor merge -i ./ -r ${catalogue} -o merged -s
    """
}

process meteor_strain {
    tag { reads_id }
    conda "meteor=2.0.21"
    memory { 
        def baseMemoryGB = params.minimum_memory.toGiga()
        def memoryGB = params.fast ? Math.min(baseMemoryGB, 10) : baseMemoryGB
        "${memoryGB}G"
    }

    input:
    tuple val(reads_id), path(mapping)
    path(catalogue)

    output:
    path("strain/*"), emit: strains, optional: true

    """
    meteor strain -i ${mapping} -r ${catalogue} -o strain
    """
}

process meteor_tree {
    cpus params.cpus
    conda "meteor=2.0.21"
    publishDir "$myDir", mode: 'copy'

    input:
    path(strain)

    output:
    path("tree")

    """
    meteor tree -i  ./ -r  -o tree -t ${params.cpus}
    """
}

workflow {
    if (params.catalogue_name) {
        catalogue_ch = meteor_download().catalogue
    } else if (params.catalogue != "") {
        catalogue_ch = Channel.value(file(params.catalogue))
    } else {
        exit 1, "ERROR: Either --catalogue_name or --catalogue must be provided"
    }

    readChannel = Channel.fromFilePairs("${params.in}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}", flat: true)
                    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}"}
                    .map { sample_id, file1, file2 ->
                        def count = file1.countFastq() / 500000
                        [sample_id, file1, file2, count.round(2)]
                    }
    meteor_fastq(readChannel)
    meteor_mapping(meteor_fastq.out, catalogue_ch)
    meteor_profile(meteor_mapping.out.mapping, catalogue_ch)
    profiles = meteor_profile.out.profile.map { id, profpath -> profpath}
    collected_prof = profiles.collect()
    meteor_merge(collected_prof, catalogue_ch)
    meteor_strain(meteor_mapping.out.mapping, catalogue_ch)
    strains = meteor_strain.out.strains.collect(flat: false)
    meteor_tree(strains)
}

    