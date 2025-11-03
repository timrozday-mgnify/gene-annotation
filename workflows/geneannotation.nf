include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SEQKIT_SEQ } from '../modules/nf-core/seqkit/seq/main'
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


    // Filter out very long sequences

    SEQKIT_SEQ(cdss)
    filtered_cdss = SEQKIT_SEQ.out.fastx


    // Annotate CDSs

    // Pfam
    // split, comine with db, group/count then flatten again, and do a final map
    if(!params.skip_pfam) {
        chunked_cdss_pfam_in = filtered_cdss
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
    }

    // dbCAN3
    // split, comine with db, group/count then flatten again, and do a final map
    if(!params.skip_dbcan) {
        chunked_cdss_dbcan_in = filtered_cdss
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
    }

    // KOfam
    // split, comine with db, group/count then flatten again, and do a final map
    if(!params.skip_kofam) {
        chunked_cdss_kofam_in = filtered_cdss
            .splitFasta(
                size: params.kofam_chunksize,
                elem: 1,
                file: true
            )

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

        kofam_seqs_dbs_ch = chunked_cdss_kofam_in
            .combine(kofam_chunked_db_ch)
            .map { meta, seqs, db_meta -> [meta, [seqs, db_meta]] }
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
            .multiMap { meta, v ->
                def (seqs, db_meta) = v
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
    }

    emit:
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}