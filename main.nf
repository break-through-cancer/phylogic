nextflow.enable.dsl=2

params.maf       = null
params.seg       = null
params.purity    = null
params.timepoint = 1
params.outdir    = "results"

process ADD_LOCAL_CN_TO_MAF {
    input:
    tuple path(maf), path(seg), val(purity), val(timepoint)

    output:
    tuple path("phylogic_input.maf"), path(seg), val(purity), val(timepoint)

    script:
    """
    set -euo pipefail

    python3 <<'PY'
import csv

maf_in = "${maf.getName()}"
seg_in = "${seg.getName()}"
maf_out = "phylogic_input.maf"

def norm_name(x):
    return x.strip().replace("\\ufeff", "").lower()

def norm_chrom(x):
    x = str(x).strip()
    return x[3:] if x.lower().startswith("chr") else x

segments = []
with open(seg_in, "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\\t")
    fields = reader.fieldnames or []
    fmap = {norm_name(c): c for c in fields}

    chrom_col = fmap.get("chromosome") or fmap.get("chrom")
    start_col = fmap.get("start.bp") or fmap.get("loc.start") or fmap.get("start")
    end_col   = fmap.get("end.bp") or fmap.get("loc.end") or fmap.get("end")

    a1_col = fmap.get("modal.a1") or fmap.get("rescaled.cn.a1") or fmap.get("expected.a1")
    a2_col = fmap.get("modal.a2") or fmap.get("rescaled.cn.a2") or fmap.get("expected.a2")

    if not all([chrom_col, start_col, end_col, a1_col, a2_col]):
        raise SystemExit("ERROR: Could not identify required SEG columns. Found columns: " + ",".join(fields))

    for row in reader:
        segments.append((
            norm_chrom(row[chrom_col]),
            int(float(row[start_col])),
            int(float(row[end_col])),
            row[a1_col],
            row[a2_col]
        ))

comments = []
header = None
rows = []

with open(maf_in, "r", newline="") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            comments.append(line)
            continue
        header = line.rstrip("\\n").split("\\t")
        break

    if header is None:
        raise SystemExit("ERROR: Could not find MAF header")

    reader = csv.DictReader(f, fieldnames=header, delimiter="\\t")
    rows = list(reader)

    MAX_CLUSTER_MUTS = 10000

    if len(rows) > MAX_CLUSTER_MUTS:
      print("Downsampling MAF rows for clustering:", len(rows), "->", MAX_CLUSTER_MUTS)
      rows = rows[:MAX_CLUSTER_MUTS]

    


maf_fmap = {norm_name(c): c for c in header}
chrom_col = maf_fmap.get("chromosome")
start_col = maf_fmap.get("start_position")

if chrom_col is None or start_col is None:
    raise SystemExit("ERROR: MAF missing Chromosome or Start_position column")

out_header = list(header)
if "local_cn_a1" not in out_header:
    out_header.append("local_cn_a1")
if "local_cn_a2" not in out_header:
    out_header.append("local_cn_a2")

matched = 0
unmatched = 0

for row in rows:
    chrom = norm_chrom(row[chrom_col])
    pos = int(float(row[start_col]))

    hit = None
    for seg_chrom, seg_start, seg_end, a1, a2 in segments:
        if chrom == seg_chrom and seg_start <= pos <= seg_end:
            hit = (a1, a2)
            break

    if hit:
        row["local_cn_a1"] = hit[0]
        row["local_cn_a2"] = hit[1]
        matched += 1
    else:
        row["local_cn_a1"] = "1"
        row["local_cn_a2"] = "1"
        unmatched += 1

with open(maf_out, "w", newline="") as out:
    for line in comments:
        out.write(line)
    writer = csv.DictWriter(out, fieldnames=out_header, delimiter="\\t", extrasaction="ignore", lineterminator="\\n")
    writer.writeheader()
    writer.writerows(rows)

print("Wrote", maf_out)
print("Matched mutations:", matched)
print("Unmatched defaulted to 1/1:", unmatched)
PY
    """
}

process MAKE_SIF {
    input:
    tuple path(maf), path(seg), val(purity), val(timepoint)

    output:
    tuple path("patient_id.txt"), path("patient.sif"), path(maf), path(seg)

    script:
    """
    set -euo pipefail

    SAMPLE_ID=\$(awk -F'\\t' '
      !/^#/ && NF > 1 {
        for (i=1; i<=NF; i++) {
          if (tolower(\$i)=="sample") sample_col=i
          if (tolower(\$i)=="tumor_sample_barcode") tumor_col=i
        }
        getline
        if (sample_col) print \$sample_col
        else if (tumor_col) print \$tumor_col
        else exit 1
        exit
      }
    ' ${maf})

    PATIENT_ID=\$(echo "\$SAMPLE_ID" | awk -F'.' '{print \$1"." \$2}')

    echo "\$PATIENT_ID" > patient_id.txt

    printf "sample_id\\tmaf_fn\\tseg_fn\\tpurity\\ttimepoint\\n" > patient.sif
    printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\$SAMPLE_ID" "${maf.getName()}" "${seg.getName()}" "${purity}" "${timepoint}" >> patient.sif

    cat patient.sif
    """
}

process RUN_CLUSTER {
    container "gcr.io/broad-getzlab-workflows/phylogicndt:v50"
    publishDir "${params.outdir}/cluster", mode: "copy"

    input:
    tuple path(patient_id_txt), path(sif), path(maf), path(seg)

    output:
    tuple path(patient_id_txt), path(sif), path("cluster_out")

    script:
    """
    set -euo pipefail

    PATIENT_ID=\$(cat ${patient_id_txt})
    export PYTHONPATH=/build/PhylogicNDT:/build/PhylogicNDT/data:\${PYTHONPATH:-}

    mkdir -p cluster_out
    cp ${sif} cluster_out/
    cp ${maf} cluster_out/
    cp ${seg} cluster_out/

    cd cluster_out

    python2 /build/PhylogicNDT/PhylogicNDT.py Cluster \
      -i "\$PATIENT_ID" \
      -sif patient.sif \
      --Pi_k_mu 10 \
      --Pi_k_r 10

    ls -lah
    """
}

process RUN_BUILDTREE {
    container "gcr.io/broad-getzlab-workflows/phylogicndt:v50"
    publishDir "${params.outdir}/buildtree", mode: "copy"

    input:
    tuple path(patient_id_txt), path(sif), path(cluster_out)

    output:
    path "buildtree_out"

    script:
    """
    set -euo pipefail

    PATIENT_ID=\$(cat ${patient_id_txt})
    export PYTHONPATH=/build/PhylogicNDT:/build/PhylogicNDT/data:\${PYTHONPATH:-}

    mkdir -p buildtree_out
    cp ${sif} buildtree_out/patient.sif
    cp -r ${cluster_out}/* buildtree_out/

    cd buildtree_out

    MUTATION_CCF=\$(find . -maxdepth 1 -type f \\( -iname "*mut*ccf*" -o -iname "*mutation*ccf*" \\) | head -n 1 | sed 's#^./##')
    CLUSTER_CCF=\$(find . -maxdepth 1 -type f -iname "*cluster*ccf*" | head -n 1 | sed 's#^./##')

    echo "Using mutation CCF: \$MUTATION_CCF"
    echo "Using cluster CCF: \$CLUSTER_CCF"

    if [ -z "\$MUTATION_CCF" ] || [ -z "\$CLUSTER_CCF" ]; then
      echo "ERROR: Could not find Cluster outputs" >&2
      ls -lah >&2
      exit 1
    fi

    python2 /build/PhylogicNDT/PhylogicNDT.py BuildTree \
        -i "\$PATIENT_ID" \
        -sif patient.sif \
        -m "\$MUTATION_CCF" \
        -c "\$CLUSTER_CCF"

    ls -lah
    """
}

workflow {
    input_ch = Channel.of([
        file(params.maf),
        file(params.seg),
        params.purity,
        params.timepoint ?: 1
    ])

    maf_with_cn_ch = ADD_LOCAL_CN_TO_MAF(input_ch)
    sif_ch = MAKE_SIF(maf_with_cn_ch)
    cluster_ch = RUN_CLUSTER(sif_ch)
    RUN_BUILDTREE(cluster_ch)
}