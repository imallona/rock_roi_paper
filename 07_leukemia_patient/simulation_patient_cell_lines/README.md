# Expectations on simulations

The simulated reads have an exps field called expb which is 1 if the expectation is to detect the reads by bwa aln, while it will be 0 if we do not expect to detect the read with bwa aln. The field exps is 1 if we expect to detect the read with STAR fusion while it will be 0 if we do not expect to detect the read.<br />
The non-fused ABL and BCR is expected to be detected with bwa aln but not STAR fusion (as it is not fused). <br />
Additionally, the non-fused BCR will always multimapped as the annotated transcripts in the gencode.v38.pc_transcripts.fa have very similar sequences. 
