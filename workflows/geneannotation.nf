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
    dbs_path_ch = FETCHDB.out.dbs

    dbs_path_ch
        .branch { meta, _fp ->
            pfam: meta.id == 'pfam'
        }
        .set { dbs }


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
    pfam_db = dbs.pfam
        .map { meta, fp ->
            file("${fp}/${meta.files.hmm}")
        }
        .first()

    // split, comine with db, group/count then flatten again, and do a final map
    chunked_cdss_pfam_in = cdss
        .splitFasta(
            size: params.hmmsearch_chunksize,
            elem: 1,
            file: true
        )
        .combine(pfam_db)
        .map { meta, reads, db -> tuple(meta, tuple(reads, db)) }
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
            return [meta, db, seqs, false, true, true] 
        }

    HMMER_HMMSEARCH(chunked_cdss_pfam_in)

    CONCATENATE(
        HMMER_HMMSEARCH.out.domain_summary
        .groupTuple()
        .map{ meta, results -> tuple(meta, "${meta.id}.domtbl.gz", results) }
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
