# Aim and warning

This set of reports aims to explore, quantify and troubleshoot several aspects of our technique:

1. The abundance of offtargets (<0.5% ontargets for `mixing` [v2](https://www.biorxiv.org/content/10.1101/2024.05.18.594120v2.full))
   1. The capture of TSO-looking reads
   2. The capture of rock-looking transcripts
   3. Other sources (GC content etc)
2. The generation of chimeric reads with misprimed ROI primers (`mixing`, `leukemia`, `pdgfra`)
   1. Quantification, description of the problem
   2. The discrimination between fusion reads and chimeric reads (`leukemia`)
   3. Association to other features (local motifs, GC content etc)

These thoughts and concerns belong to Izaskun Mallona (izaskun.mallona@mls.uzh.ch) only.

Our experimental design offers the necessary flexibility to explore these aspects. Namely:
- `mixing`. Normally processed, e.g. evaluating _only_ unique mappers or raw reads, provides rock, roi, and unmod roi (in WTA), helping 1.1, 1.2, 1.3, 2.2
- `pdgfra`. Validation set for `mixing` with a single reference genome.
- `leukemia`. Cell lines have associated truths (fusions to be detected) and one of the ROIs (major) is resistant against chimeras due to its placement away from the breakpoint.

We have preliminary data describing aspects of these two major concerns but have not formalized them into a single report.

# Design

- Hypothesis-driven
- Topdown analysis on downsampled data and self-contained
- From GEO data, whenever possible - not relying on other analysis

# 2.1 Quantification of ROI potential artifacts

## Molecular biology

- `rock` allows capturing mRNAs by sequence similarity. It acts as a replacement of oligodTs (which capture any polyA transcript).
  - We hardly see the `rock` captured sequences because they typically are very 3' of the sequenced read; for short reads, we do not reach it.
- `roi` is as primer replacing random hexamers. Sequencing reads start with the ROI sequence in their 5'. `roi` primers show high sequence similarity to some region of the targets, typically 5' of the `rock` capture (if present) or of the polyA (if `rock` is not present)
  - We use very short `roi` primers and low temperature; hence, `roi` priming is expected to be somewhat unspecific
  - If that is true, a sizeable portion of the TSO and WTA reads should contain the `roi` sequences
  - If that is true, nonspecific priming events can result in `roi`-`cdna` chimeric reads where the `roi` sequence is not canonically expected to be 5' of the (unspecific) cDNA sequence
  - If that is true, specific priming events will result in `roi`-`cdna` reads where the `roi` sequence is the oligo sequence (it is not cDNA, but a primer), even if the genome/transcriptome has variants on it
  - To check this, we can make use of the unmodified experiments as well as `rockroi` and `wta roi` (no rock) experiments.
  - To check this, we can make use of the leukemia experiments where read lengths are long enough as to sequence both `rock` and `roi` parts within single reads
    - So we can build expected read conformations (= true read fusions vs chimeric reads; akin to simulations)
      - And, if these read conformations are different, we could quantify them separately (goes to aim 2.2) 

## Read conformations

Schemas depicting expected ontarget and offtarget read conformations, per experiment and perceived challenges:

## Quantification

First, we run simple ROI scans on raw reads as follows:
