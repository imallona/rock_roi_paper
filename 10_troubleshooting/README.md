# Aim

This set of reports aims to explore, quantify and troubleshoot several aspects of our technique:

1. The abundance of offtargets (<0.5% ontargets for `mixing` [v2](https://www.biorxiv.org/content/10.1101/2024.05.18.594120v2.full))
   1.1 The capture of TSO-looking reads
   1.2 The capture of rock-looking transcripts
   1.3 Other sources (GC content etc)
2. The generation of chimeric reads with misprimed ROI primers (`mixing`, `leukemia`, `pdgfra`)
   2.1. Quantification, description of the problem
   2.2. The discrimination between fusion reads and chimeric reads (`leukemia`)
   2.2. Association to other features (local motifs, GC content etc)


Our experimental design offers the necessary flexibility to explore these aspects. Namely:
- `mixing`. Normally processed, e.g. evaluating _only_ unique mappers or raw reads, provides rock, roi, and unmod roi (in WTA), helping 1.1, 1.2, 1.3, 2.2
- `pdgfra`. Validation set for `mixing` with a single reference genome.
- `leukemia`. Cell lines have associated truths (fusions to be detected) and one of the ROIs (major) is resistant against chimeras due to its placement away from the breakpoint.

# Design

- Hypothesis-driven
- Topdown analysis on downsampled data and self-contained
- From GEO data, whenever possible

# 1.1

# 1.2 

# 1.3

# 2.1 Quantification of ROI potential artifacts

First, we run simple ROI scans on raw reads as follows:
