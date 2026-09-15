#!/usr/bin/env nextflow

params.help=false
params.in = ""
params.out = ""
params.cpus =  4
params.fast = false
params.check_catalogue = false
params.catalogue_name = ""
params.catalogue = ""
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

process meteor_download {
    tag { params.catalogue_name }
    conda "meteor=2.0.22"
    
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
    conda "meteor=2.0.22"

    input:
    tuple val(reads_id), path(forward), path(reverse), val(count)

    output:
    tuple val(reads_id), path("fastq/*"), val(count)

    script:
    """
    meteor fastq -i ./ -p -o fastq
    """
}

process meteor_mapping {
    tag { reads_id }
    cpus params.cpus
    memory {
        // Base memory scales with the catalogue size: complete catalogues need 3x
        // their size in GB, light (taxo) catalogues only 1/3. One extra GB per
        // million reads is added on top.
        def base_mem = 3 * cat_size
        if (cat_type == 'taxo') base_mem = cat_size / 3
        def reads_mem = Math.ceil((count as int) / 1_000_000)
        return (base_mem + reads_mem).GB
    }
    conda "meteor=2.0.22"

    input:
    tuple val(reads_id), path(fastq), val(count)
    tuple path(catalogue), val(cat_size), val(cat_type)

    output:
    tuple val(reads_id), path("mapping/*"), emit: mapping

    script:
    """
    meteor mapping -i ${fastq} -r ${catalogue} -t ${params.cpus} -o mapping --kf
    """
}

process meteor_profile {
    tag { reads_id }
    cpus params.cpus
    memory {
        // Profiling loads the catalogue: 1x its size for complete catalogues,
        // 1/3 for light (taxo) catalogues.
        def base_mem = cat_size
        if (cat_type == 'taxo') base_mem = cat_size / 3
        return base_mem.GB
    }
    conda "meteor=2.0.22"

    input:
    tuple val(reads_id), path(mapping)
    tuple path(catalogue), val(cat_size), val(cat_type)

    output:
    tuple val(reads_id), path("profile/*"), emit: profile

    script:
    """
    meteor profile -i ${mapping} -r ${catalogue} -o profile
    """
}

process meteor_merge {
    memory {
        // Memory grows with the number of samples profiled: 20MB per sample per
        // catalogue GB (2MB for taxo catalogues), plus a 1/10 catalogue size floor.
        def sample_count = (profile instanceof List ? profile : [profile])
                .collect { p -> p.parent.toString() }.unique().size()
        def slope = (cat_type == 'taxo') ? (cat_size / 500) : (cat_size / 50)
        def intercept = cat_size / 10
        return (slope * sample_count + intercept).GB
    }
    conda "meteor=2.0.22"
    publishDir params.out, mode: 'copy'

    input:
    path(profile)
    tuple path(catalogue), val(cat_size), val(cat_type)

    output:
    path("merged")

    script:
    """
    meteor merge -i ./ -r ${catalogue} -o merged -s
    """
}

process meteor_strain {
    tag { reads_id }
    conda "meteor=2.0.22"
    memory {
        // Strain profiling is the most demanding step: 4x catalogue size plus a
        // 5 GB floor for complete catalogues, 0.4x for taxo catalogues.
        def mem = 5 + 4 * cat_size
        if (cat_type == 'taxo') mem = 5 + 0.4 * cat_size
        return mem.GB
    }

    input:
    tuple val(reads_id), path(mapping)
    tuple path(catalogue), val(cat_size), val(cat_type)

    output:
    path("strain/*"), emit: strains, optional: true

    script:
    """
    meteor strain -i ${mapping} -r ${catalogue} -o strain
    """
}

process meteor_tree {
    cpus params.cpus
    conda "meteor=2.0.22"
    memory {
        // Tree inference scales with the number of samples analysed.
        def sample_count = (strain instanceof List ? strain : [strain])
                .collect { p -> p.parent.toString() }.unique().size()
        return (0.2 * sample_count + 20).GB
    }
    publishDir params.out, mode: 'copy'

    input:
    path(strain)

    output:
    path("tree")

    script:
    """
    meteor tree -i  ./ -r  -o tree -t ${params.cpus}
    """
}


workflow {

    // Parameter validation and setup
    def allowed_catalogues = params.allowed_catalogues.split(',')

    if(params.help){
        usage()
        exit(1)
    }

    // Validate catalogue_name parameter
    if (params.catalogue_name && !allowed_catalogues.contains(params.catalogue_name)) {
        println "ERROR: Invalid catalogue name '${params.catalogue_name}'"
        println "Allowed catalogues are:"
        allowed_catalogues.each { catalogue -> println "  - ${catalogue}" }
        exit 1
    }

    file(params.out).mkdirs()

    // Reference catalogue sizing, used to scale process memory allocations:
    //   - catalogue_size is parsed from the catalogue name (e.g. hs_10_4_gut -> 10.4)
    //   - database_type is 'taxo' in fast mode (light catalogue), otherwise 'complete'
    //   - with a custom --catalogue, both are read from the reference.json file
    def catalogue_size = null
    def database_type = 'complete'

    if (params.catalogue_name) {
        def matcher = (params.catalogue_name =~ /^[a-z]+_([0-9]+)_([0-9]+)_/)
        if (!matcher) {
            println "ERROR: Unexpected catalogue name format '${params.catalogue_name}'"
            exit 1
        }
        catalogue_size = "${matcher[0][1]}.${matcher[0][2]}".toDouble()
        database_type = params.fast ? 'taxo' : 'complete'
        catalogue_ch = meteor_download().catalogue
    } else if (params.catalogue != "") {
        def json_file = file(params.catalogue).listFiles()?.find { f -> f.name.endsWith('reference.json') }
        if (!json_file) {
            println "ERROR: No *_reference.json found in ${params.catalogue}"
            exit 1
        }
        def meta = new groovy.json.JsonSlurper().parseText(json_file.text)
        def ref_name = meta.reference_info.reference_name
        database_type = meta.reference_info.database_type ?: 'complete'
        def matcher = (ref_name =~ /^[a-z]+_([0-9]+)_([0-9]+)_/)
        if (!matcher) {
            println "ERROR: Unexpected reference name format '${ref_name}'"
            exit 1
        }
        catalogue_size = "${matcher[0][1]}.${matcher[0][2]}".toDouble()
        catalogue_ch = channel.value(file(params.catalogue))
    } else {
        exit 1, "ERROR: Either --catalogue_name or --catalogue must be provided"
    }

    // Attach catalogue size and database type so process memory directives can use them
    catalogue_ch = catalogue_ch.map { c -> tuple(c, catalogue_size, database_type) }
    
        readChannel = channel.fromFilePairs("${params.in}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}", flat: true)
                        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.in}"}
                        .map { sample_id, file1, file2 ->
                            def count = file1.countFastq()
                            [sample_id, file1, file2, count]
                        }
        meteor_fastq(readChannel)
        meteor_mapping(meteor_fastq.out, catalogue_ch)
        meteor_profile(meteor_mapping.out.mapping, catalogue_ch)
        profiles = meteor_profile.out.profile.map { _id, profpath -> profpath }
        collected_prof = profiles.collect()
        meteor_merge(collected_prof, catalogue_ch)
        meteor_strain(meteor_mapping.out.mapping, catalogue_ch)
        strains = meteor_strain.out.strains.collect(flat: false)
        meteor_tree(strains)

    workflow.onComplete = {
        // any workflow property can be used here
        println "Pipeline complete"
        println "Command line: $workflow.commandLine"
    }

    workflow.onError = {
        println "Oops .. something went wrong"
    }
}
