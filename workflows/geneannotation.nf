/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { HMMER_HMMSEARCH } from '../modules/nf-core/hmmer/hmmsearch/main'
include { RUNDBCAN_CAZYMEANNOTATION } from '../modules/nf-core/rundbcan/cazymeannotation/main'
include { KOFAMSCAN } from '../modules/nf-core/kofamscan/main'
include { CONCATENATE as PFAM_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as DBCAN_OVERVIEW_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as DBCAN_CAZYME_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as DBCAN_SUBSTRATE_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as DBCAN_DIAMOND_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as KOFAM_TSV_CONCATENATE } from  '../modules/local/concatenate/main'
include { CONCATENATE as KOFAM_TXT_CONCATENATE } from  '../modules/local/concatenate/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENEANNOTATION {
    main:
    ch_versions = Channel.empty()

    // Parse samplesheet and fetch reads
    samplesheet = Channel.fromList(samplesheetToList(params.samplesheet, "${workflow.projectDir}/assets/schema_input.json"))

    cdss = samplesheet.map {
        id, faa ->
        [
            ['id': id],
            file(faa),
        ]
    }

    // Annotate CDSs

    // Pfam
    // split, comine with db, group/count then flatten again, and do a final map
    chunked_cdss_pfam_in = cdss
        .splitFasta(
            size: params.pfam_chunksize,
            elem: 1,
            file: true
        )
        .groupTuple()
        .flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            def chunksize = chunks_.size()
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunksize), chunk)
            }
        }
        .map { meta, seqs ->
            return [
                meta, 
                file(params.databases.pfam.files.profiles), 
                params.databases.pfam.n_profiles, 
                seqs, 
                false, true, true
            ] 
        }

    HMMER_HMMSEARCH(chunked_cdss_pfam_in)

    PFAM_CONCATENATE(
        HMMER_HMMSEARCH.out.domain_summary
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_${meta.db_id}.domtbl.gz", results) }
    )

    // dbCAN3
    // split, comine with db, group/count then flatten again, and do a final map
    chunked_cdss_dbcan_in = cdss
        .splitFasta(
            size: params.dbcan_chunksize,
            elem: 1,
            file: true
        )
        .groupTuple()
        .flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            def chunksize = chunks_.size()
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunksize), chunk)
            }
        }

    RUNDBCAN_CAZYMEANNOTATION(chunked_cdss_dbcan_in, file(params.databases.dbcan.files.db))

    DBCAN_OVERVIEW_CONCATENATE(
        RUNDBCAN_CAZYMEANNOTATION.out.cazyme_annotation
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_dbcan_overview.tsv", results) }
    )
    DBCAN_CAZYME_CONCATENATE(
        RUNDBCAN_CAZYMEANNOTATION.out.dbcanhmm_results
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_dbcan_cazyme.tsv", results) }
    )
    DBCAN_SUBSTRATE_CONCATENATE(
        RUNDBCAN_CAZYMEANNOTATION.out.dbcansub_results
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_dbcan_substrate.tsv", results) }
    )
    DBCAN_DIAMOND_CONCATENATE(
        RUNDBCAN_CAZYMEANNOTATION.out.dbcandiamond_results
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_dbcan_diamond.tsv", results) }
    )

    // KOfam
    // split, comine with db, group/count then flatten again, and do a final map
    chunked_cdss_kofam_in = cdss
        .splitFasta(
            size: params.kofam_chunksize,
            elem: 1,
            file: true
        )
    chunked_cdss_kofam_in.view { "chunked_cdss_kofam_in - ${it}" }

    kofam_chunked_db_ch = Channel
        .from(
            params.databases.kofam_chunked.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('files')) {
                        return [id: k] + v
                    }
                }
            }
        )
    kofam_chunked_db_ch.view { "kofam_chunked_db_ch - ${it}" }

    kofam_seqs_dbs_ch = chunked_cdss_kofam_in
        .combine(kofam_chunked_db_ch)
        .groupTuple()
        .flatMap {
            meta, chunks ->
            def chunks_ = chunks instanceof Collection ? chunks : [chunks]
            def chunksize = chunks_.size()
            return chunks_.collect {
                chunk ->
                tuple(groupKey(meta, chunksize), chunk)
            }
        }
        .multiMap { meta, seqs, db_meta ->
            seqs: [meta, seqs]
            profiles: file(db_meta.files.profiles)
            ko_list: file(db_meta.files.ko_list)
        }
    
    KOFAMSCAN(
        kofam_seqs_dbs_ch.seqs,
        kofam_seqs_dbs_ch.profiles,
        kofam_seqs_dbs_ch.ko_list,
    )
    KOFAM_TSV_CONCATENATE(
        KOFAMSCAN.out.tsv
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_kofam.tsv", results) }
    )
    KOFAM_TXT_CONCATENATE(
        KOFAMSCAN.out.txt
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_kofam.txt", results) }
    )
    emit:
    cds_locations = cdss
    pfam_annotations = PFAM_CONCATENATE.out.concatenated_file
    dbcan_overview_annotations = DBCAN_OVERVIEW_CONCATENATE.out.concatenated_file
    dbcan_cazyme_annotations = DBCAN_CAZYME_CONCATENATE.out.concatenated_file
    dbcan_substrate_annotations = DBCAN_SUBSTRATE_CONCATENATE.out.concatenated_file
    dbcan_diamond_annotations = DBCAN_DIAMOND_CONCATENATE.out.concatenated_file
    kofam_tsv_annotations = KOFAM_TSV_CONCATENATE.out.concatenated_file
    kofam_txt_annotations = KOFAM_TXT_CONCATENATE.out.concatenated_file
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
