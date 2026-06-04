nextflow.enable.dsl=2

params.samplesheet = null

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
include { SAMTOOLS_FAIDX }   from './modules/nf-core/samtools/main'


process FILTER_VARIANTS {

    tag "$sample_id"

    publishDir "${params.outdir}/filtered_variants",
               mode: 'copy'

    input:
    tuple val(group), val(sample_id), path(vcf)

    output:
    tuple val(group), val(sample_id),
          path("${sample_id}.filtered.vcf")

    script:
    """
    bcftools filter \
        -i 'QUAL>20' \
        $vcf \
        -Ov \
        -o ${sample_id}.filtered.vcf
    """

    stub:
    """
    touch ${sample_id}.filtered.vcf
    """
}


workflow {

    if (!params.samplesheet) {
        error "Provide --samplesheet"
    }

    /*
     * Read all samples into one channel
     */

    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->

            tuple(
                row.group,
                row.sample_id,
                file(row.r1),
                file(row.r2),
                file(row.reference)
            )
        }

    /*
     * Split by group field
     */

    grouped_samples_ch = samples_ch.groupTuple(by: 0)

    /*
     * Join back to one channel
     */

    reads_ch = grouped_samples_ch
        .flatMap { group, samples ->

            samples.collect { s ->

                tuple(
                    group,
                    s[1],
                    s[2],
                    s[3],
                    s[4]
                )
            }
        }

    /*
     * FASTQC input format:
     * tuple(sample_id, r1, r2)
     */

    fastqc_input = reads_ch.map {
        group, sample_id, r1, r2, ref ->

        tuple(sample_id, r1, r2)
    }

    FASTQC_RAW(fastqc_input)

    trimmed_ch = TRIM(fastqc_input)

    FASTQC_TRIMMED(trimmed_ch)

    /*
     * Recover group information after trimming
     */

    trim_info_ch = reads_ch.map {
        group, sample_id, r1, r2, ref ->

        tuple(sample_id, group, ref)
    }

    trimmed_with_group = trimmed_ch
        .join(trim_info_ch)
        .map { sample_id, r1, r2, group, ref ->

            tuple(
                group,
                sample_id,
                r1,
                r2,
                ref
            )
        }

    /*
     * Mapping
     */

    map_input = trimmed_with_group.map {
        group, sample_id, r1, r2, ref ->

        tuple(sample_id, r1, r2)
    }

    ref_ch = reads_ch
        .map { group, sample_id, r1, r2, ref -> ref }
        .unique()

    mapped_ch = MAP(map_input, ref_ch)

    PLOT_COVERAGE(mapped_ch)

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
                reference,
                fai
            )
        }

    /*
     * Add group information back to BAM files
     */

    sample_group_ch = reads_ch.map {
        group, sample_id, r1, r2, ref ->

        tuple(sample_id, group)
    }

    bam_for_variants = mapped_ch
        .join(sample_group_ch)
        .map { sample_id, bam, bai, group ->

            tuple(
                [id: sample_id],
                bam,
                bai
            )
        }

    /*
     * Variant calling
     */

    mpileup_ch = BCFTOOLS_MPILEUP(
        bam_for_variants,
        reference_with_index,
        false
    )

    call_ch = BCFTOOLS_CALL(
        mpileup_ch.out,
        false
    )

    /*
     * Join back to one channel
     */

    variants_ch = call_ch.out.vcf
        .join(sample_group_ch)
        .map { sample_id, vcf, group ->

            tuple(
                group,
                sample_id,
                vcf
            )
        }

    /*
     * Final analysis
     */

    FILTER_VARIANTS(variants_ch)
}
