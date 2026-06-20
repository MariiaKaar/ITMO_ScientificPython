nextflow.enable.dsl=2

// =====================
// PARAMS
// =====================
params.reads     = null
params.reference = null

// =====================
// FASTQC
// =====================
process FASTQC {

    tag "$sample_id ($stage)"
    cpus 8

    input:
    val stage
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    fastqc -t ${task.cpus} $r1 $r2
    """
}

// =====================
// TRIM
// =====================
process TRIM {

    tag "$sample_id"
    cpus 8

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id),
          path("trimmed_1.fastq"),
          path("trimmed_2.fastq")

    script:
    """
    trimmomatic PE -threads ${task.cpus} \
        $r1 $r2 \
        trimmed_1.fastq unpaired_1.fastq \
        trimmed_2.fastq unpaired_2.fastq \
        SLIDINGWINDOW:4:20 MINLEN:50
    """
}

// =====================
// ASSEMBLY (fallback)
// =====================
process ASSEMBLE {

    tag "assembly"
    cpus 8

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "contigs.fasta"

    script:
    """
    spades.py -1 $r1 -2 $r2 -o spades_out -t ${task.cpus}
    cp spades_out/contigs.fasta contigs.fasta
    """
}

// =====================
// INDEX REFERENCE
// =====================
process INDEX_REF {

    tag "index"
    cpus 8

    input:
    path ref

    output:
    path "ref.*"

    script:
    """
    bwa index $ref
    """
}

// =====================
// MAPPING
// =====================
process MAP {

    tag "$sample_id"
    cpus 8

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref

    output:
    tuple val(sample_id), path("aligned.bam")

    script:
    """
    bwa mem -t ${task.cpus} $ref $r1 $r2 | samtools sort -@ ${task.cpus} -o aligned.bam
    """
}

// =====================
// COVERAGE + PLOT
// =====================
process PLOT_COVERAGE {

    tag "$sample_id"
    cpus 8

    conda 'bioconda::samtools=1.19 conda-forge::python=3.10 conda-forge::matplotlib=3.8.0'
    publishDir "results/coverage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id),
          path("${sample_id}_coverage.png"),
          path("${sample_id}_depth.txt")

    script:
    """
    samtools depth $bam > ${sample_id}_depth.txt

    cat <<EOF > plot.py
import matplotlib.pyplot as plt

depths = []
with open("${sample_id}_depth.txt") as f:
    for line in f:
        depths.append(int(line.strip().split()[2]))

plt.figure(figsize=(12,4))
plt.plot(depths, linewidth=0.5)
plt.title("Coverage: ${sample_id}")
plt.xlabel("Position")
plt.ylabel("Depth")
plt.grid(True)
plt.tight_layout()
plt.savefig("${sample_id}_coverage.png", dpi=300)
EOF

    python3 plot.py
    """
}

// =====================
// FASTQC WRAPPERS (fix DSL2 reuse bug)
// =====================
workflow FASTQC_RAW {
    take:
    reads_ch

    main:
    FASTQC("raw", reads_ch)
}

workflow FASTQC_TRIMMED {
    take:
    reads_ch

    main:
    FASTQC("trimmed", reads_ch)
}

// =====================
// MAIN WORKFLOW (NAMED)
// =====================
workflow main_pipeline {

    // -----------------
    // INPUT (LOCAL FASTQ ONLY)
    // -----------------
    if (!params.reads) {
        error "Provide --reads (paired FASTQ files)"
    }

    reads_ch = Channel.fromPath(params.reads)

    // -----------------
    // QC RAW
    // -----------------
    FASTQC_RAW(reads_ch)

    // -----------------
    // TRIM
    // -----------------
    trimmed_ch = TRIM(reads_ch)

    // -----------------
    // QC TRIMMED
    // -----------------
    FASTQC_TRIMMED(trimmed_ch)

    // -----------------
    // REFERENCE OR ASSEMBLY
    // -----------------
    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference)
    } else {
        ref_ch = ASSEMBLE(trimmed_ch)
    }

    // -----------------
    // INDEX
    // -----------------
    INDEX_REF(ref_ch)

    // -----------------
    // MAP
    // -----------------
    mapped = MAP(trimmed_ch, ref_ch)

    // -----------------
    // COVERAGE PLOT
    // -----------------
    PLOT_COVERAGE(mapped)
}