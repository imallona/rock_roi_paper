# Aim and disclaimer

This set of reports aims to explore, quantify and troubleshoot several aspects of our technique:

1. The abundance of offtargets (<0.5% ontargets for `mixing` [v2](https://www.biorxiv.org/content/10.1101/2024.05.18.594120v2.full))
   1. The capture of transcripts with `tso` motifs
   2. The capture of transcripts with `rock` motifs
   3. Other sources (GC content etc)
2. The generation of chimeric reads with misprimed `roi` primers (`mixing`, `leukemia`, `pdgfra`)
   1. Quantification, description of the problem
   2. The discrimination between fusion reads and chimeric reads (`leukemia`)
   3. Association to other features (local motifs, GC content etc)
   4. Impact in cancer applications, e.g. BCR/ABL profiling
   
These thoughts are written by Izaskun Mallona (izaskun.mallona@mls.uzh.ch) **only**. Our research is meant to be reproducible and robust; this report might go too deep into further exploring possible caveats when processing clinical samples in a very specific setting (gene fusions). This is probably a consequence of being a cancer researcher.

Our experimental design offers the necessary flexibility to explore these features and possible caveats. Namely:
- `mixing`. Normally processed, e.g. evaluating _only_ unique mappers or raw reads, provides rock, roi, and unmod roi (in WTA) experiments, helping 1.1, 1.2, 1.3, 2.2
- `pdgfra`. Validation set for `mixing` with a single reference genome; low capture efficiency if I recall correctly
- `leukemia`. Cell lines have associated known truths (fusions to be detected) and one of the roi primers (major) is easy to tell apart from chimeras due to its placement away from the breakpoint. Not so for other roi primers.

We have preliminary data describing aspects of these two major concerns but have not formalized them into a single report nor manuscript section.

# Design

We plan our checks to be:

- Molecular biology-driven (hypothesis-driven). Not hypothesis-free nor exploratory. Wet lab protocols and data properties are linked, and we can intuitively expect some challenges and data properties from the former.
- Topdown analysis on downsampled data and self-contained, so not linked to other PRs (e.g. leukemia, nor leukemia_patient)
- From GEO data, whenever possible - not relying on other analysis.

# 2.1 Quantification of ROI potential artifacts

## Naming

rock, ROCK, RoCK etc mean both technique and capture sequence; same for ROI, roi etc. (roi|rock)`m` highlights sequence _motifs_ are being discussed, e.g. sequences present within reads. 

## Molecular biology

- `rock` allows capturing mRNAs by sequence similarity. It acts as a replacement of oligodTs (which capture polyA transcript).
  - We hardly see the `rockm` (rock motifs) within reads because they typically downstream of it; for short reads, we do not reach it. We have longer sequencing results for `leukemia` though.
- `roi` is as primer replacing random hexamers. Sequencing reads start with the `roi` sequence in their 5'. `roi` primers show high sequence similarity to some region of the targets. They are 5' of the `rock` capture (if present) or of the polyA (if `rock` is not present)
  - We use *very short* `roi` primers and low annealing temperature; hence, `roi` priming is expected to be unspecific. This is a feature, not a bug. But it generates offtargets, more precisely chimeric reads.
  - If true, a sizeable portion of the TSO and WTA reads should contain `roi` sequences at the read start, if running `roi`. Not if running `rock` alone.
  - If true, nonspecific priming events can result in `roi`-`cdna` chimeric reads where the `roi` sequence was not expected to be 5' of the (unspecific) cDNA sequence
  - If true, specific priming events will result in `roi`-`cdna` reads where the `roi` sequence is the oligo sequence (it is not cDNA, but a primer), even if the genome/transcriptome has variants on it
  - To check this:
    - We can make use of the unmodified experiments as well as `rockroi` and `wta roi` (no rock) experiments.
    - We can make use of the leukemia experiments where read lengths are long enough as to sequence both `rock` and `roi` parts within single reads
       - So we can build expected read conformations (= true read fusions vs chimeric reads; akin to simulations)
          - And, if these read conformations are different, we can quantify them separately (goes to aim 2.2) 

## Read conformations

Schemas depicting expected ontarget and offtarget read conformations, per experiment and perceived challenges:

- "|" means nucleotide identity
- "*" means mismatch
- `R` depicts a ROIm nucleotide. Mind these are very short!
- `C` depicts a RoCKm nucleotide. Mind these are longer.

### Random primed polyA

Sequencing errors aside, cDNA read and mRNA are identical

```
-------------------------------AAAAA mRNA (not seen)
||||||||||||||||||||||               identity
----------------------               cDNA read
```

### ROI primed polyA ontarget

3' of the read is cDNA; 5' is the ROI primer. Requirements: distance between ROIm (ROI motif) and polyA is short; seq variation in ROIs (ROI sequence) is only technical.

```
RRR----------------------------AAAAA mRNA (not seen)
||||||||||||||||||||||               identity
RRR-------------------               cDNA read
```

### ROI primed polyA chimeric

The 5' of the read is congruent with a ROIs but the mRNA does not have a ROIm; the read is hence chimeric and alignable with softclipping.

```
-------------------------------AAAAA mRNA (not seen)
***|||||||||||||||||||               identity
RRR------------------$               cDNA read
```

### Random primed rock

Unless sequenced long enought the ROCKm (rock motif) is not present in the read.

```
----------------------------CCCCC----- mRNA (not seen)
||||||||||||||||||||||                 identity
----------------------                 cDNA read
```

### ROI primed rock ontarget

Unless sequenced long enough the ROCKm (rock motif) is not present in the read.

```
RRR-------------------------CCCCC----- mRNA (not seen)
||||||||||||||||||||||                 identity
RRR-------------------                 cDNA read
```

### ROI primed rock offtarget

Not sure why is it captured; due to the TSO? pure WTA abundance? Worth checking distance to polyA? 

```
RRR----------------------------------- mRNA (not seen)
||||||||||||||||||||||                 identity
RRR-------------------                 cDNA read
```

### ROI primed rock offtarget chimeric

Neither the ROCKm nor ROIm are present in the mRNA. Alignable with softclipping aligners.

```
-------------------------------------- mRNA (not seen)
***|||||||||||||||||||                 identity
RRR-------------------                 cDNA read
```

### Challenge: fused transcripts vs chimeric reads

Lets imagine a gene `B` and gene `A` which are fused into `BA`.

```
--------------------------------------   B mRNA (not seen)
        ||||||||||||******************   B vs BA identity
        ------------------------------   BA fused mRNA (not seen)
              ******||||||||||||||||||   BA vs A identity
              ------------------------   A mRNA (not seen)
```

This can be approached by designing `rock` rather downstream of `A` and `roi` upstream of the breakpoint and in `B`. Hence capturing `B` unfused and `BA` fused, but guiding the reads to the fused `BA` (via `roi`).

```
--------------------------------CCCC--   B mRNA (not seen)
        ||||||||||||******************   B vs BA identity
        ------------------------------   BA fused mRNA (not seen)
              ******||||||||||||||||||   BA vs A identity
              RRRR--------------------   A mRNA (not seen)
```

So the read conformation is :

```
        <----B-----><----A----------->   unfused B and A origins
               -->  <--                  tolerance
        ------RRRR--------------CCCC--   BA fused mRNA (not seen)
              |||||||||||||||            identity
              RRRR-----------            cDNA read
```

Problem is, the tolerance region between the 3' of the roim and the fusion point length **matters**: this region allows to distinguish fusions from ROI chimeras.

Tolerance calculation: sequence long enough so the `rockm` capture is also part of the read, so true cDNAs and chimeras can be told apart:

- True fusions contain the expected BA slice starting with `roim` and ending with (part of) `rockm`
- Chimeras are shaped as `rockm` : slice of `A` of length (read_length - roi_length) upstream of the rockm start position. With potential variants, including insertions and deletions. Chimeric sequence and insertions and deletions depend on the (unfused) `A` isoforms expressed and (perhaps) their expression level.

```
        <----B-----><----A----------->   unfused B and A origins
               -->  <--                  tolerance
        ------RRRR--------------CCCC--   BA fused mRNA (not seen)
              ||||||||||||||||||||       identity
              RRRR--------------CC       cDNA read (true fusion)
              RRRR**||||||||||||CC       expected chimera
```

An alternative view, BA fused transcripts could be:
- oligodT, randomly primed. Unstacked reads, perhaps showing the fusion (unlikely, 3' bias)
- oligodT, roi. Stacked reads, unlikely due to distance between oligodT and roi
- rock, roi. Stacked reads, read starts at roi seq and goes over the breakpoint and `A` until read length
- chimeras

To explore these, diagnostic aspects are: roim presence; variants within the tolerance region; abundances as compared to other background artifacts, e.g.

- ROI:B quantities
- ROI:anything chimeras

As extra datasets and validations, we will reanalyze:

- Truths where the ROI primer is >20 nt away from the fusion point

## Quantification

First, we run simple ROI scans on raw reads as follows:


### Truths, experimental metadata

Cell lines:

|   Cell line       | Fusion   | ROIs  |
| :---------------- | :------: | ----: |
| a                 |   b      | c     |

# Analysis

## `extract_and_quantify_roi_reads` ROIm quantifier

Script with signature `extract_and_quantify_roi_reads --roi regex --cdna cdna.fastq.gz --cb cbumi.fastq.gz --output annotated.fastq.gz`

1. Decorates cDNA and fastqs with UMI and CB (from R1) and whether WTA or TSO (by adding these information to the read id)
2. Quantifies and writes fastq records containing the roi pattern(s) to the output file

Applications:
- To all experiments, downsampled. Why? To report the ROI rates.

Modifications:
- Perhaps add a 'trim' function to remove the ROI part to facilitate downstream analysis

## `align_roi_reads`

Aligns and counts `extract_and_quantify_roi_reads` outputs, perhaps with the `trim` function, to the relevant genome/transcriptome. Signature `align_roi_reads --fastq annotated.fastq.gz --genome x`.

Applications:
- First crack at quantifying offtarget/chimeras.
- All experiments, downsampled

