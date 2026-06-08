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

include { SAMTOOLS_FAIDX }    from './modules/nf-core/samtools/main'

process FILTER_VARIANTS {

    tag "$sample_id"

    publishDir "results/filtered_variants", mode: 'copy'

    input:
    tuple val(group), val(sample_id), path(vcf)

    output:
    tuple val(group),
          val(sample_id),
          path("${sample_id}.filtered.vcf.gz")

    script:
    """
    bcftools filter \
        -i 'QUAL>20' \
        $vcf \
        -Oz \
        -o ${sample_id}.filtered.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.filtered.vcf.gz
    """
}

workflow {

    if (!params.samplesheet) {
        error "Provide --samplesheet"
    }

    if (!params.reference) {
        error "Provide --reference"
    }

    /*
     * Read all samples into one channel
     */

    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row ->
            tuple(
                row.group,
                row.sample_id,
                file(row.r1),
                file(row.r2)
            )
        }

    /*
     * Split by group
     */

    grouped_ch = samples_ch.groupTuple(by:0)

    /*
     * Join back to one channel
     */

    reads_ch = samples_ch

    /*
     * Remove group for existing processes
     */

    reads_for_pipeline = reads_ch.map {
        group, sample_id, r1, r2 ->
        tuple(sample_id, r1, r2)
    }

    /*
     * Save group information
     */

    sample_group_ch = samples_ch.map {
        group, sample_id, r1, r2 ->
        tuple(sample_id, group)
    }

    ref_ch = Channel.fromPath(params.reference)

    FASTQC_RAW(reads_for_pipeline)

    trimmed_ch = TRIM(reads_for_pipeline)

    FASTQC_TRIMMED(trimmed_ch)

    mapped_ch = MAP(trimmed_ch, ref_ch)
    mapped_ch.view()

    PLOT_COVERAGE(
        mapped_ch.map { sample_id, bam, bai ->
            tuple(sample_id, bam)
        }
    )

    /*
     * Reference indexing
     */

    reference_for_index = ref_ch.map { ref ->

        tuple(
            [id: 'reference'],
            ref,
            []
        )

    }

    SAMTOOLS_FAIDX(reference_for_index, false)

    reference_with_index = ref_ch
        .join(SAMTOOLS_FAIDX.out.fai)
        .map { reference, fai ->
            tuple(
            [id: 'reference'], reference, fai)
        }

    /*
     * Prepare BAMs for mpileup
     */

    bam_for_variants = mapped_ch.map { sample_id, bam, bai ->

        tuple(
            [id: sample_id],
            bam,
            [],
            []
        )

    }
    bam_for_variants.view()

    /*
     * Variant calling
     */

    mpileup_result = BCFTOOLS_MPILEUP(
        bam_for_variants,
        reference_with_index,
        false
    )
    mpileup_result.vcf.view()

    variants_ch = mpileup_result.vcf
    .combine(sample_group_ch)
    .map { vcf_tuple, group_tuple ->
        def (meta, vcf) = vcf_tuple
        def (sample_id, group) = group_tuple
        tuple(group, sample_id, vcf)
    }
    variants_ch.view {"VARIANTS: $it"}


    FILTER_VARIANTS(variants_ch)

}
