nextflow.enable.dsl=2

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

process MAP {

    tag "$sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(r1), path(r2)
    path ref

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai")

    script:
    """
    bwa index $ref

    bwa mem -t ${task.cpus} $ref $r1 $r2 | \
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam

    samtools index ${sample_id}.sorted.bam
    """
}

process PLOT_COVERAGE {

    tag "$sample_id"

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

x = []
with open("${sample_id}_depth.txt") as f:
    for line in f:
        x.append(int(line.strip().split()[2]))

plt.figure(figsize=(12,4))
plt.plot(x, linewidth=0.5)
plt.title("Coverage: ${sample_id}")
plt.xlabel("Position")
plt.ylabel("Depth")
plt.tight_layout()
plt.savefig("${sample_id}_coverage.png", dpi=300)
EOF

    python plot.py
    """
}

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
