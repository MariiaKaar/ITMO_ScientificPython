nextflow.enable.dsl=2

params.reads     = null
params.reference = null

// =====================
// FASTQC
// =====================
process FASTQC {

    tag "$sample_id ($stage)"
    cpus 4

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
    cpus 4

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
// MAP (FIXED)
// =====================
process MAP {

    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    bwa index $ref
    bwa mem -t ${task.cpus} $ref $r1 $r2 | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    """
}

// =====================
// COVERAGE + PLOT
// =====================
process PLOT_COVERAGE {

    tag "$sample_id"
    cpus 4

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
plt.tight_layout()
plt.savefig("${sample_id}_coverage.png", dpi=300)
EOF

    python3 plot.py
    """
}

// =====================
// WORKFLOWS
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
// MAIN
// =====================
workflow main_pipeline {

    if (!params.reads) {
        error "Provide --reads"
    }

    reads_ch = Channel
        .fromFilePairs(params.reads,checkIfExists: true)
        .map { sample_id, reads ->
            tuple(sample_id, reads[0], reads[1])}

    FASTQC_RAW(reads_ch)

    trimmed_ch = TRIM(reads_ch)

    FASTQC_TRIMMED(trimmed_ch)

    if (!params.reference) {
        error "Provide --reference (.fna supported)"
    }

    ref_ch = Channel.fromPath(params.reference)



    mapped = MAP(trimmed_ch, ref_ch)

    PLOT_COVERAGE(mapped)
}
