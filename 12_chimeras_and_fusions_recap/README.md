When quantifying gene fusions, we recommend not placing the ROI too close to the breakpoint because ROIseq is prone to producing (mispriming) chimeric reads and teasing apart true chimeras from technical chimeras can be difficult.

Our [debugging process is open](https://github.com/imallona/rock_roi_paper/pull/20), including [issue 11](https://github.com/imallona/rock_roi_method/issues/11) (for fusions) and [some results](https://github.com/imallona/rock_roi_paper/tree/main/10_troubleshooting).

*If aiming to quantify gene fusions please mind [point 66 of the wet lab protocol](https://www.protocols.io/view/rock-and-roi-bead-modification-library-generation-rm7vzjyb5lx1/v2?step=66)* to design ROI primers in a way biological chimeras (gene fusions) are clearly distinct from mispriming chimeras. Also see the [general README](https://github.com/imallona/rock_roi_method?tab=readme-ov-file#tldr-sequencing-reads-layout) to get an idea of read layout from ROI-seq experiments.

We have a (perhaps) private repository on debugging chimeric reads at [rock_roi_debug](https://github.com/imallona/rock_roi_debug).

Our paper describes a BCR::ABL1 detection on a very suboptimal experimental design. [Its analysis](https://github.com/imallona/rock_roi_debug/blob/main/08_agreed_flow.sh) is not recommended not necessary if you place the ROI primers further apart from the breakpoint. Our analysis includes:

1. extract reads with a roi primer, not anchored, and trim the left side, while keeping the roi sequence. Hence, these reads must: have the breakpoint to be counted later, and be long enough as to be mappable
2. decorate reads with the valid/allowedlisted cellbarcodes (allowlist from running the standard pipeline on WTA using an unfused normal genome/transcriptome, e.g. Gencode's); discard reads from cells not matching the CB allowedlist
3. align reads from `2` to a citable fusion db, such as [ChimerDB 4.0](https://www.kobic.re.kr/chimerdb/download)
4. count reads (alignments) overlapping the fusion and without deletions/softclipping

We recommend using `cutadapt` for steps 1 and 2.

We suggest using `STAR` with disabled softclipping and splicing for step 3, e.g.

```bash
## align with STAR all reads fulfilling the two filtering steps above: that is, ROI-containing, and having
##   a relaxedly valid WTA/TSO barcode
## do not allow softclipping (endtoend alignment)
STAR --runThreadN $NTHREADS \
     --genomeDir "$WD"/../06_align_to_a_fusion_database_with_star/fusions_no_splicing_idx/ \
     --outTmpDir "tmp_star" \
     --readFilesCommand zcat \
     --outFileNamePrefix "star_fusion/" \
     --readFilesIn labelled_cdna.fq.gz \
     --quantMode TranscriptomeSAM GeneCounts \
     --sjdbGTFfile ../cherrypicked_fusions.gtf \
     --twopassMode Basic \
     --alignEndsType EndToEnd \
     --outSAMunmapped Within \
     --outSAMtype BAM Unsorted \
     --alignIntronMax 1 \
     --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 &> star_fusion.log
```

We recommend grafting the CB as readgroups and counting reads *that overlap the fusion breakpoint* using `featureCounts` or similar.
