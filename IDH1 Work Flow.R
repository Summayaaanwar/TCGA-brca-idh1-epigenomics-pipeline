# ============================================================
#  workflow figure (Wet-lab clinical cohort + TCGA-BRCA validation)
# - Exports: SVG (vector), PDF (vector), PNG (high-res)
# - Larger fonts + spacing to avoid overlaps
# ============================================================

# ---- Packages ----
if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
  install.packages("DiagrammeR", dependencies = TRUE)
}
library(DiagrammeR)

# For export (recommended for journals)
if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) {
  install.packages("DiagrammeRsvg", dependencies = TRUE)
}
if (!requireNamespace("rsvg", quietly = TRUE)) {
  install.packages("rsvg", dependencies = TRUE)
}

library(DiagrammeRsvg)
library(rsvg)

# ---- Controls for publication export ----
out_prefix   <- "IDH1_workflow"  # output file prefix
target_width_in <- 7.5           # typical single-column ~3.5in, double-column ~7–7.5in
dpi_equiv <- 600                 # common journal raster requirement
png_width_px <- as.integer(target_width_in * dpi_equiv)

# ---- DOT (Graphviz) ----
workflow_dot <- "
digraph workflow {

  // -------------------------
  // Global styling (publication)
  // -------------------------
  graph [
    rankdir=TB,
    layout=dot,
    splines=ortho,
    newrank=true,

    // more spacing for readability at larger fonts
    nodesep=0.50,
    ranksep=0.85,
    pad=0.18,
    margin=0.10,

    // larger global font
    fontsize=18,
    fontname='Helvetica',

    // helps compact without crowding
    ratio=compress
  ];

  node [
    shape=box,
    style='rounded,filled',
    fillcolor='white',
    color='gray35',
    penwidth=1.25,

    // larger node font
    fontsize=15,
    fontname='Helvetica',

    // a bit more padding so text doesn't touch borders
    margin='0.22,0.14'
  ];

  edge [
    color='gray35',
    penwidth=1.10,
    arrowsize=0.85
  ];

  // -------------------------
  // Clinical cohort (top)
  // -------------------------
  A [
    label='Clinical cohort\\nBreast tumor + matched ANCT\\nEthics approval + informed consent\\nClinical variables recorded',
    penwidth=1.45
  ];

  // =========================================================
  // Wet-lab (clinical cohort)
  // =========================================================
  subgraph cluster_wetlab {
    label='Wet-lab (clinical cohort)';
    labelloc='t';
    labeljust='l';

    // larger cluster label
    fontsize=16;
    fontname='Helvetica';

    style='rounded,dashed';
    color='gray55';
    penwidth=1.20;

    B  [label='Tissue processing\\nAliquot for DNA and RNA'];

    // DNA arm
    C1 [label='DNA extraction'];
    C2 [label='Bisulfite conversion\\n(converted DNA)'];
    C3 [label='MSP (M/U primers)\\nTarget regions: I1, I2, distal element'];
    C4 [label='Agarose gel electrophoresis\\n+ imaging (fixed settings)'];
    C5 [label='ImageJ densitometry\\nBand intensity: M and U'];
    C6 [label='Derived methylation indices\\nPM, M-value, ΔMeth\\n(semi-quantitative)'];

    // RNA arm
    R1 [label='RNA extraction'];
    R2 [label='cDNA synthesis'];
    R3 [label='qRT-PCR\\nIDH1 expression (normalized)\\n+ QC/replicates'];

    // Wet-lab flow
    A  -> B;
    B  -> C1;
    C1 -> C2 -> C3 -> C4 -> C5 -> C6;

    B  -> R1 -> R2 -> R3;
  }

  // =========================================================
  // Statistical analysis (clinical cohort)
  // =========================================================
  S [
    label='Statistical analysis (clinical cohort)\\nPaired tests, group comparisons\\nAssociations with clinicopathologic variables\\nSurvival (Kaplan–Meier/Cox where applicable)',
    penwidth=1.45
  ];

  C6 -> S;
  R3 -> S;

  // =========================================================
  // In silico validation (TCGA-BRCA)
  // =========================================================
  subgraph cluster_tcga {
    label='In silico validation (TCGA-BRCA)';
    labelloc='t';
    labeljust='l';

    fontsize=16;
    fontname='Helvetica';

    style='rounded,dashed';
    color='gray55';
    penwidth=1.20;

    T1 [label='TCGA-BRCA data retrieval\\n450K methylation + RNA-seq\\nClinical + survival annotation'];
    T2 [label='Preprocessing + QC\\nfiltering/normalization'];

    T3c [label='Stratification\\n(IDH1-high/low)'];
    T3a [label='DMP analysis\\n(450K methylation)'];
    T3b [label='DGE analysis\\n(RNA-seq)'];
    T3d [label='Survival models\\n(Kaplan–Meier/Cox on TCGA)'];

    // TCGA flow
    T1 -> T2 -> T3c;

    // Stratification feeds all analyses
    T3c -> T3a;
    T3c -> T3b;
    T3c -> T3d;

    // Align analysis boxes (cleaner)
    { rank=same; T3a; T3b; T3d; }
  }

  // =========================================================
  // Integration
  // =========================================================
  I [
    label='Integration\\nCompare clinical cohort with TCGA\\nConsistency checks + interpretation',
    penwidth=1.70
  ];

  S   -> I;
  T3a -> I;
  T3b -> I;
  T3d -> I;

  // =========================================================
  // Genome-wide context (AFTER Integration as requested)
  // =========================================================
  G [
    label='Genome-wide context\\n100 kb probe-density clustering\\nAnnotation + enrichment\\nGene–CpG network/hubs',
    penwidth=1.45
  ];
  I -> G;

  // =========================================================
  // Outputs
  // =========================================================
  O [
    label='Outputs\\nMain figures + Supplementary (gels, ImageJ)\\nReproducible code + README\\nData/Code availability statements',
    penwidth=1.45
  ];
  G -> O;

  // Align merge contributors row
  { rank=same; S; T3a; T3b; T3d; }
}
"

# ---- Render ----
viz <- grViz(workflow_dot)
viz

# ---- Export: SVG (vector) ----
svg_txt <- export_svg(viz)
svg_file <- paste0(out_prefix, ".svg")
writeLines(svg_txt, con = svg_file)

# ---- Export: PDF (vector via rsvg) ----
pdf_file <- paste0(out_prefix, ".pdf")
rsvg_pdf(svg_file, file = pdf_file)

# ---- Export: PNG (high-res fallback) ----
png_file <- paste0(out_prefix, "_", dpi_equiv, "dpi.png")
rsvg_png(svg_file, file = png_file, width = png_width_px)

