/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { samplesheetToList } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { SEQSTATS } from  '../modules/local/seqstats/main'
include { CONCATENATE } from  '../modules/local/concatenate/main'
include { PYRODIGAL as PYRODIGAL_SMALL } from '../modules/nf-core/pyrodigal/main'
include { PYRODIGAL as PYRODIGAL_LARGE } from '../modules/nf-core/pyrodigal/main'
include { HMMER_HMMSEARCH } from '../modules/nf-core/hmmer/hmmsearch/main'
include { FETCHDB } from '../subworkflows/local/fetchdb/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENEANNOTATION {
    main:
    ch_versions = Channel.empty()

    // Fetch databases
    db_ch = Channel
        .from(
            params.databases.collect { k, v ->
                if ((v instanceof Map) && v.containsKey('files')) {
                    return [id: k] + v
                }
            }
        )
        .filter { it }

    FETCHDB(db_ch, "${projectDir}/${params.databases.cache_path}")
    hmm_dbs = FETCHDB.out.dbs
        .map { meta, fp ->
        def db_params = params.databases[meta.id] 
        [meta, [file("${fp}/${db_params.files.hmm}"), db_params.variables.num_models]] }

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

    // split, comine with db, group/count then flatten again, and do a final map
    chunked_cdss_pfam_in = cdss
        .splitFasta(
            size: params.hmmsearch_chunksize,
            elem: 1,
            file: true
        )
        .combine(hmm_dbs)
        .map { 
            meta, reads, db_meta, db -> tuple(meta + ['db_id': db_meta.id], tuple(reads, db)) }
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
        .map { meta, v ->
            def (seqs, db) = v
            def (db_fp, db_nseqs) = db
            return [meta, db_fp, db_nseqs, seqs, false, true, true] 
        }

    HMMER_HMMSEARCH(chunked_cdss_pfam_in)

    CONCATENATE(
        HMMER_HMMSEARCH.out.domain_summary
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}_${meta.db_id}.domtbl.gz", results) }
    )

    emit:
    cds_locations = cdss
    functional_annotations = CONCATENATE.out.concatenated_file
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
