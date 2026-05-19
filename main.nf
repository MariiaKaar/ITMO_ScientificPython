nextflow.enable.dsl=2

params.reads     = null
params.reference = null

include {
    FASTQC
    FASTQC_RAW
    FASTQC_TRIMMED
    TRIM
    MAP
    PLOT_COVERAGE
} from './modules/local/processes'

include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_CALL }    from './modules/nf-core/bcftools/call/main'
include { SAMTOOLS_FAIDX }    from './modules/nf-core/samtools/main'

workflow {

    if (!params.reads) {
        error "Provide --reads"
    }

    if (!params.reference) {
        error "Provide --reference"
    }

    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads ->
            tuple(sample_id, reads[0], reads[1])
        }

    ref_ch = Channel.fromPath(params.reference)


    FASTQC_RAW(reads_ch)

    trimmed_ch = TRIM(reads_ch)

    FASTQC_TRIMMED(trimmed_ch)

    mapped_ch = MAP(trimmed_ch, ref_ch)


    PLOT_COVERAGE(mapped_ch)

    reference_for_index = ref_ch.map { ref ->

        tuple([id: 'reference'], ref, [])

    }


    SAMTOOLS_FAIDX(reference_for_index, false)


    reference_with_index = ref_ch

        .join(SAMTOOLS_FAIDX.out.fai)

        .map { reference, fai ->

            tuple(reference, fai)

        }


    bam_for_variants = mapped_ch.map { sample_id, bam, bai ->

        tuple(

            [id: sample_id],

            bam,

            bai

        )

    }


//    BCFTOOLS_MPILEUP(

//        bam_for_variants,

//        reference_with_index,

//        false

//    )

}
