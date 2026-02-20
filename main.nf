nextflow.enable.dsl = 2

def required = ['store_root', 'project_name', 'input_plink', 'bed_file', 'reference_genome']

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
    tag "chr${chr}"

    publishDir "${params.store_root}/${params.project_name}/temp/step_03_slice",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chr),
          path(pgen, stageAs: "in.pgen"),
          path(pvar, stageAs: "in.pvar"),
          path(psam, stageAs: "in.psam"),
          path(bed_file)

    output:
    tuple val(chr),
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
    publishDir "${params.store_root}/${params.project_name}/temp/step_04_merge",
        mode: 'copy',
        overwrite: true

    input:
    path sliced_files

    output:
    path "step04_merge_list.txt", emit: merge_list
    path "step04_merged.pgen", emit: pgen
    path "step04_merged.pvar", emit: pvar
    path "step04_merged.psam", emit: psam

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
    publishDir "${params.store_root}/${params.project_name}/temp/step_05_detect_duplicates",
        mode: 'copy',
        overwrite: true

    input:
    path merged_pvar
    path dup_detector_script

    output:
    path "step05_exclude_variants.txt", emit: exclude
    path "step05_detect_duplicates.log", emit: log

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
    publishDir "${params.store_root}/${params.project_name}/temp/step_06_clean",
        mode: 'copy',
        overwrite: true

    input:
    path merged_pgen
    path merged_pvar
    path merged_psam
    path exclude_list

    output:
    path "step06_cleaned.pgen", emit: pgen
    path "step06_cleaned.pvar", emit: pvar
    path "step06_cleaned.psam", emit: psam
    path "step06_cleaned.log", emit: log

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
    publishDir "${params.store_root}/${params.project_name}/temp/step_07_normalize",
        mode: 'copy',
        overwrite: true

    input:
    path step6_pgen
    path step6_pvar
    path step6_psam
    path reference_fasta
    path reference_fasta_fai

    output:
    path "step07_normalized.vcf"
    tuple path("step07_final.pgen"), path("step07_final.pvar"), path("step07_final.psam"), emit: normalized
    path "step07_final.log"

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
    tag "${score_name}"

    publishDir "${params.store_root}/${params.project_name}/temp/step_08_score",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(score_name), path(score_file), path(norm_pgen), path(norm_pvar), path(norm_psam)

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

    def sliceResults = Channel
        .from(params.input_plink)                  // params.input_plink is a List
        .map { it.toString().trim() }
        .filter { it }                              // drop blanks just in case
        .map { entry ->
            // Accept "prefix" or "*.pgen/*.pvar/*.psam" (strip extension if present)
            def prefix = entry.replaceFirst(/\.(pgen|pvar|psam)$/, '')

            // Extract chromosome label from the prefix (e.g. chr20 -> "20")
            def m = (prefix =~ /.*chr([0-9XYM]+).*/)
            def chrLabel = m.matches() ? m[0][1] : 'NA'
            def safeChrLabel = chrLabel.replaceAll(/[^A-Za-z0-9._-]/, '_')

            tuple(
                safeChrLabel,
                file("${prefix}.pgen", checkIfExists: true),
                file("${prefix}.pvar", checkIfExists: true),
                file("${prefix}.psam", checkIfExists: true),
                bed
            )
        }
        | STEP03_SLICE_CHR

    def step4 = sliceResults
        .filter { chr, outPgen, outPvar, outPsam -> outPgen && outPvar && outPsam }
        .map    { chr, outPgen, outPvar, outPsam -> [outPgen, outPvar, outPsam] }
        .flatten()
        .collect()
        .filter { files -> files && files.size() > 0 }
        | STEP04_MERGE_SLICES

    def step5 = STEP05_DETECT_DUPLICATES(step4.pvar, dupDetectorScript)

    def step6 = STEP06_CLEAN(step4.pgen, step4.pvar, step4.psam, step5.exclude)

    def step7 = STEP07_NORMALIZE(step6.pgen, step6.pvar, step6.psam, hg38Ref, hg38RefFai)

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

        scoreFiles
            .combine(step7.normalized)
            .map { scoreName, scoreFile, normPgen, normPvar, normPsam -> tuple(scoreName, scoreFile, normPgen, normPvar, normPsam) }
            | STEP08_SCORE
    } else {
        log.info "No score_file provided; skipping STEP08_SCORE."
    }
}
