nextflow.enable.dsl = 2

def required = ['store_root', 'project_name', 'ancestries', 'input_plink_pattern', 'bed_file', 'reference_genome']

def isBlank = { v ->
    if (v == null) return true
    if (v instanceof Collection)
        return v.isEmpty() || v.every { (it == null) || it.toString().trim() == '' }
    return v.toString().trim() == ''
}

def missing = required.findAll { key -> isBlank(params[key]) }

if (missing) {
    error "Missing required params: ${missing.join(', ')}"
}

process STEP03_SLICE_CHR {
    tag "${ancestry}:chr${chr}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_03_slice",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry),
          val(chr),
          path(pgen, stageAs: "in.pgen"),
          path(pvar, stageAs: "in.pvar"),
          path(psam, stageAs: "in.psam"),
          path(bed_file)

    output:
    tuple val(ancestry),
      val(chr),
      path("chr${chr}.step03.slice.pgen"),
      path("chr${chr}.step03.slice.pvar"),
      path("chr${chr}.step03.slice.psam"),
      optional: true

    script:
    """
    set +e
    plink2 \\
      --pfile in \\
      --extract bed1 "${bed_file}" \\
      --var-filter  \\
      --rm-dup force-first \\
      --make-pgen \\
      --out chr${chr}.step03.slice
    rc=\$?
    set -e

    if [[ \$rc -eq 0 ]]; then
      exit 0
    fi

    if [[ \$rc -eq 13 ]] && grep -q 'No variants remaining after main filters' "chr${chr}.step03.slice.log"; then
      echo "No variants after extract for chr${chr}; skipping merge input."
      exit 0
    fi

    exit \$rc
    """
}


process STEP04_MERGE_SLICES {
    tag "${ancestry}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_04_merge",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry), path(sliced_files)

    output:
    tuple val(ancestry), path("step04_merge_list.txt"), emit: merge_list
    tuple val(ancestry), path("step04_merged.pgen"), emit: pgen
    tuple val(ancestry), path("step04_merged.pvar"), emit: pvar
    tuple val(ancestry), path("step04_merged.psam"), emit: psam

    script:
    """
    set -euo pipefail
    shopt -s nullglob

    # Collect only sliced .pgen files from staged inputs.
    pgens=( *.step03.slice.pgen )
    if (( \${#pgens[@]} == 0 )); then
      echo "No .pgen files staged for merge" >&2
      exit 1
    fi

    printf "%s\\n" "\${pgens[@]}" | sed 's/\\.pgen\$//' | sort -V > step04_merge_list.txt

    plink2 \\
      --pmerge-list step04_merge_list.txt \\
      --make-pgen \\
      --out step04_merged
    """
}


process STEP05_DETECT_DUPLICATES {
    tag "${ancestry}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_05_detect_duplicates",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry), path(merged_pvar)
    path dup_detector_script

    output:
    tuple val(ancestry), path("step05_exclude_variants.txt"), emit: exclude
    tuple val(ancestry), path("step05_detect_duplicates.log"), emit: log

    script:
    """
    set -euo pipefail
    python3 "${dup_detector_script}" --pvar "${merged_pvar}" > step05_detect_duplicates.log 2>&1
    if [[ -f step04_merged_exclude.txt ]]; then
      mv step04_merged_exclude.txt step05_exclude_variants.txt
    elif [[ -f merged.step4_exclude.txt ]]; then
      mv merged.step4_exclude.txt step05_exclude_variants.txt
    else
      echo "Duplicate detector did not create expected exclude file." >&2
      exit 1
    fi
    """
}


process STEP06_CLEAN {
    tag "${ancestry}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_06_clean",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry), path(merged_pgen), path(merged_pvar), path(merged_psam), path(exclude_list)

    output:
    tuple val(ancestry), path("step06_cleaned.pgen"), emit: pgen
    tuple val(ancestry), path("step06_cleaned.pvar"), emit: pvar
    tuple val(ancestry), path("step06_cleaned.psam"), emit: psam
    tuple val(ancestry), path("step06_cleaned.log"), emit: log

    script:
    """
    set -euo pipefail
    plink2 \\
      --pfile step04_merged \\
      --exclude "${exclude_list}" \\
      --make-pgen \\
      --out step06_cleaned
    """
}


process STEP07_NORMALIZE {
    tag "${ancestry}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_07_normalize",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry), path(step6_pgen), path(step6_pvar), path(step6_psam)
    path reference_fasta
    path reference_fasta_fai

    output:
    tuple val(ancestry), path("step07_normalized.vcf")
    tuple val(ancestry), path("step07_final.pgen"), path("step07_final.pvar"), path("step07_final.psam"), emit: normalized
    tuple val(ancestry), path("step07_final.log")

    script:
    """
    set -euo pipefail

    plink2 \\
      --pfile step06_cleaned \\
      --export vcf id-delim='^' \\
      --out step07_from_plink

    vcf_for_norm="step07_from_plink.vcf"
    if awk 'NR==1 { exit (\$1 ~ /^chr/) ? 0 : 1 }' "${reference_fasta_fai}"; then
      cat > step07_chr_rename.txt <<'EOF'
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
M chrM
EOF
      bcftools annotate \\
        --rename-chrs step07_chr_rename.txt \\
        -Ov \\
        -o step07_from_plink.renamed.vcf \\
        step07_from_plink.vcf
      vcf_for_norm="step07_from_plink.renamed.vcf"
    fi

    bcftools norm \\
      -f "${reference_fasta}" \\
      -Ov \\
      -o step07_normalized.vcf \\
      "\${vcf_for_norm}"

    plink2 \\
      --vcf step07_normalized.vcf \\
      --id-delim '^' \\
      --set-all-var-ids 'chr@:#:\$r:\$a' \\
      --make-pgen \\
      --out step07_final
    """
}


process STEP08_SCORE {
    tag "${ancestry}:${score_name}"

    publishDir "${params.store_root}/${params.project_name}/${ancestry}/temp/step_08_score",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(ancestry), val(score_name), path(score_file), path(norm_pgen), path(norm_pvar), path(norm_psam)

    output:
    path "${score_name}.sscore"
    path "${score_name}.raw"
    path "${score_name}.score.log"
    path "${score_name}.export.log"
    path "${score_name}.sscore.vars"

    script:
    """
    set -euo pipefail

    plink2 \\
      --pfile step07_final \\
      --score "${score_file}" 1 2 3 list-variants \\
      --out "${score_name}"
    cp -f "${score_name}.log" "${score_name}.score.log"

    plink2 \\
      --pfile step07_final \\
      --export A \\
      --export-allele "${score_file}" \\
      --extract "${score_file}" \\
      --out "${score_name}"
    cp -f "${score_name}.log" "${score_name}.export.log"
    """
}


workflow {
    def bed = file(params.bed_file, checkIfExists: true)
    def hg38Ref = file(params.reference_genome, checkIfExists: true)
    def hg38RefFai = file("${params.reference_genome}.fai", checkIfExists: true)
    def dupDetectorScript = file("${projectDir}/code/detect_dup_typed_imputed.py", checkIfExists: true)

    // Fan out: one entry per (ancestry, chr)
    def sliceResults = Channel
        .from(params.ancestries)
        .flatMap { anc ->
            (1..22).collect { chr ->
                def prefix = params.input_plink_pattern
                    .replace('{ancestry}', anc)
                    .replace('{chr}', chr.toString())
                tuple(
                    anc,
                    chr.toString(),
                    file("${prefix}.pgen", checkIfExists: true),
                    file("${prefix}.pvar", checkIfExists: true),
                    file("${prefix}.psam", checkIfExists: true),
                    bed
                )
            }
        }
        | STEP03_SLICE_CHR

    // Group sliced files by ancestry, then merge per ancestry
    def step4 = sliceResults
        .filter { anc, chr, outPgen, outPvar, outPsam -> outPgen && outPvar && outPsam }
        .map    { anc, chr, outPgen, outPvar, outPsam -> tuple(anc, [outPgen, outPvar, outPsam]) }
        .groupTuple()
        .map    { anc, fileLists -> tuple(anc, fileLists.flatten()) }
        | STEP04_MERGE_SLICES

    def step5 = STEP05_DETECT_DUPLICATES(step4.pvar, dupDetectorScript)

    // Join pgen/pvar/psam/exclude by ancestry for STEP06
    def step6input = step4.pgen
        .join(step4.pvar)
        .join(step4.psam)
        .join(step5.exclude)
        .map { anc, pgen, pvar, psam, excl -> tuple(anc, pgen, pvar, psam, excl) }
    def step6 = STEP06_CLEAN(step6input)

    // Join cleaned pgen/pvar/psam by ancestry for STEP07
    def step7input = step6.pgen
        .join(step6.pvar)
        .join(step6.psam)
        .map { anc, pgen, pvar, psam -> tuple(anc, pgen, pvar, psam) }
    def step7 = STEP07_NORMALIZE(step7input, hg38Ref, hg38RefFai)

    def scoreEntries = []
    if (!isBlank(params.score_file)) {
        if (params.score_file instanceof Collection) {
            scoreEntries = params.score_file
                .collectMany { it.toString().split(',') as List }
                .collect { it.trim() }
                .findAll { it }
        } else {
            scoreEntries = params.score_file
                .toString()
                .split(',')
                .collect { it.trim() }
                .findAll { it }
        }
    }

    if (scoreEntries) {
        def scoreFiles = Channel
            .from(scoreEntries.unique())
            .map { p ->
                def base = p.tokenize('/').last()
                def scoreName = base.replaceFirst(/\\.score$/, '')
                tuple(scoreName, file(p, checkIfExists: true))
            }

        // Combine each ancestry's normalized plink with all score files
        step7.normalized
            .combine(scoreFiles)
            .map { anc, normPgen, normPvar, normPsam, scoreName, scoreFile ->
                tuple(anc, scoreName, scoreFile, normPgen, normPvar, normPsam)
            }
            | STEP08_SCORE
    } else {
        log.info "No score_file provided; skipping STEP08_SCORE."
    }
}
