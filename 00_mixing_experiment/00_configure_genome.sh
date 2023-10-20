#!/bin/bash
##
## Generates a mouse + human + gfp + mouse genome to analyze the mixinge experiment
##
## Started 19th Oct 2023
## Izaskun Mallona


WD=~/ebrunner_spectral/mixing
mkdir -p $WD/genomes
cd $WD/genomes


tab="$(printf '\t')"

cat <<EOF > captured.gtf
egfp${tab}captured${tab}exon${tab}1${tab}2818${tab}.${tab}+${tab}.${tab}gene_id "regulatory_region_egfp"; transcript_id "regulatory_region_egfp";
egfp${tab}captured${tab}exon${tab}2819${tab}2871${tab}.${tab}+${tab}.${tab}gene_id "5_utr_egfp_wpre"; transcript_id "5_utr_egfp_wpre";
egfp${tab}captured${tab}exon${tab}3593${tab}4744${tab}.${tab}+${tab}.${tab}gene_id "3_utr_egfp_wpre"; transcript_id "3_utr_egfp_wpre";
egfp${tab}captured${tab}exon${tab}2872${tab}3183${tab}.${tab}+${tab}.${tab}gene_id "5_to_roi_egfp"; transcript_id "5_to_roi_egfp";
egfp${tab}captured${tab}exon${tab}3184${tab}3245${tab}.${tab}+${tab}.${tab}gene_id "roi_egfp"; transcript_id "roi_egfp";
egfp${tab}captured${tab}exon${tab}3246${tab}3566${tab}.${tab}+${tab}.${tab}gene_id "rock_egfp"; transcript_id "rock_egfp";
egfp${tab}captured${tab}exon${tab}3567${tab}3593${tab}.${tab}+${tab}.${tab}gene_id "capture_sequence_double_egfp_egfp"; transcript_id "capture_sequence_double_egfp_egfp";
tdtomato${tab}captured${tab}exon${tab}1${tab}2818${tab}.${tab}+${tab}.${tab}gene_id "regulatory_region_tdtomato"; transcript_id "regulatory_region_tdtomato";
tdtomato${tab}captured${tab}exon${tab}2819${tab}2859${tab}.${tab}+${tab}.${tab}gene_id "5_utr_tdtomato_wpre"; transcript_id "5_utr_tdtomato_wpre";
tdtomato${tab}captured${tab}exon${tab}4291${tab}5427${tab}.${tab}+${tab}.${tab}gene_id "3_utr_tdtomato_wpre"; transcript_id "3_utr_tdtomato_wpre";
tdtomato${tab}captured${tab}exon${tab}2860${tab}3129${tab}.${tab}+${tab}.${tab}gene_id "5_to_roi_tdtomato"; transcript_id "5_to_roi_tdtomato";
tdtomato${tab}captured${tab}exon${tab}3130${tab}3192${tab}.${tab}+${tab}.${tab}gene_id "roi_1_tdtomato"; transcript_id "roi_1_tdtomato";
tdtomato${tab}captured${tab}exon${tab}3283${tab}3345${tab}.${tab}+${tab}.${tab}gene_id "roi_2_tdtomato"; transcript_id "roi_2_tdtomato";
tdtomato${tab}captured${tab}exon${tab}3856${tab}3918${tab}.${tab}+${tab}.${tab}gene_id "roi_3_tdtomato"; transcript_id "roi_3_tdtomato";
tdtomato${tab}captured${tab}exon${tab}4009${tab}4071${tab}.${tab}+${tab}.${tab}gene_id "roi_4_tdtomato"; transcript_id "roi_4_tdtomato";
tdtomato${tab}captured${tab}exon${tab}3193${tab}3282${tab}.${tab}+${tab}.${tab}gene_id "between_roi_1"; transcript_id "between_roi_1";
tdtomato${tab}captured${tab}exon${tab}3346${tab}3855${tab}.${tab}+${tab}.${tab}gene_id "between_roi_2"; transcript_id "between_roi_2";
tdtomato${tab}captured${tab}exon${tab}3919${tab}4008${tab}.${tab}+${tab}.${tab}gene_id "between_roi_3"; transcript_id "between_roi_3";
tdtomato${tab}captured${tab}exon${tab}4072${tab}4265${tab}.${tab}+${tab}.${tab}gene_id "rock_tdtomato"; transcript_id "rock_tdtomato";
tdtomato${tab}captured${tab}exon${tab}4266${tab}4290${tab}.${tab}+${tab}.${tab}gene_id "capture_sequence_double_egfp_tdtom"; transcript_id "capture_sequence_double_egfp_tdtom";
EOF

cat <<EOF > captured.fa
>egfp
GTGGCGCCCGAACAGGGACTTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGC
GCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGG
TGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAA
ATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAG
AAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACA
GTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCA
AAACAAAAGTAAGACCACCGCACAGCAAGCGGCCGCTGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGA
AGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCA
GAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGT
CAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAG
GCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCT
AAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTT
GGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGC
TTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGC
AAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAG
GTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCAC
CTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCG
ATTAGTGAACGGATCGGCACTGCGTGCGCCAATTCTGCAGACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGG
GGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACA
AATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCAGTTTGGTTAATTAACCCGTGTCGGCTCC
AGATCTGGCCTCCGCGCCGGGTTTTGGCGCCTCCCGCGGGCGCCCCCCTCCTCACGGCGAGCGCTGCCACGTCAGACGAA
GGGCGCAGCGAGCGTCCTGATCCTTCCGCCCGGACGCTCAGGACAGCGGCCCGCTGCTCATAAGACTCGGCCTTAGAACC
CCAGTATCAGCAGAAGGACATTTTAGGACGGGACTTGGGTGACTCTAGGGCACTGGTTTTCTTTCCAGAGAGCGGAACAG
GCGAGGAAAAGTAGTCCCTTCTCGGCGATTCTGCGGAGGGATCTCCGTGGGGCGGTGAACGCCGATGATTATATAAGGAC
GCGCCGGGTGTGGCACAGCTAGTTCCGTCGCAGCCGGGATTTGGGTCGCGGTTCTTGTTTGTGGATCGCTGTGATCGTCA
CTTGGTGAGTAGCGGGCTGCTGGGCTGGCCGGGGCTTTCGTGGCCGCCGGGCCGCTCGGTGGGACGGAAGCGTGTGGAGA
GACCGCCAAGGGCTGTAGTCTGGGTCCGCGAGCAAGGTTGCCCTGAACTGGGGGTTGGGGGGAGCGCAGCAAAATGGCGG
CTGTTCCCGAGTCTTGAATGGAAGACGCTTGTGAGGCGGGCTGTGAGGTCGTTGAAACAAGGTGGGGGGCATGGTGGGCG
GCAAGAACCCAAGGTCTTGAGGCCTTCGCTAATGCGGGAAAGCTCTTATTCGGGTGAGATGGGCTGGGGCACCATCTGGG
GACCCTGACGTGAAGTTTGTCACTGACTGGAGAACTCGGTTTGTCGTCTGTTGCGGGGGCGGCAGTTATGGCGGTGCCGT
TGGGCAGTGCACCCGTACCTTTGGGAGCGCGCGCCCTCGTCGTGTCGTGACGTCACCCGTTCTGTTGGCTTATAATGCAG
GGTGGGGCCACCTGCCGGTAGGTGTGCGGTAGGCTTTTCTCCGTCGCAGGACGCAGGGTTCGGGCCTAGGGTAGGCTCTC
CTGAATCGACAGGCGCCGGACCTCTGGTGAGGGGAGGGATAAGTGAGGCGTCAGTTTCTTTGGTCGGTTTTATGTACCTA
TCTTCTTAAGTAGCTGAAGCTCCGGTTTTGAACTATGCGCTCGGGGTTGGCGAGTGTGTTTTGTGAAGTTTTTTAGGCAC
CTTTTGAAATGTAATCATTTGGGTCAATATGTAATTTTCAGTGTTAGACTAGTAAATTGTCCGCTAAATTCTGGCCGTTT
TTGGCTTTTTTGTTAGACGAAGCTTGGGCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGGTCGCCACCATGGTGAGC
AAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGT
GTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGC
CCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGAC
TTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCG
CGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACA
TCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAG
GTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGG
CGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCG
ATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAAGCGGCCGC
GACTCTAGAATTCGATATCAAGCTTATCGATAATCAACCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTA
ACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTC
ATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGT
GTGCACTGTGTTTGCTGACGCAACCCCCACTGGTTGGGGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTT
TCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGCCTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACT
GACAATTCCGTGGTGTTGTCGGGGAAATCATCGTCCTTTCCTTGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGG
GACGTCCTTCTGCTACGTCCCTTCGGCCCTCAATCCAGCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTC
TTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCATCGATACCGTCGACCTCGA
GACCTAGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGATTGTGCCTGGCTAGAAGCACAAGA
GGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCC
ACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTAC
CACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTG
CTACAAGCTAGTACCAGTTGAGCAAGAGAAGGTAGAAGAAGCCAATGAAGGAGAGAACACCCGCTTGTTACACCCTGTGA
GCCTGCATGGGATGGATGACCCGGAGAGAGAAGTATTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCC
CGAGAGCTGCATCCGGACTGTACT
>tdtomato
GTGGCGCCCGAACAGGGACTTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGC
GCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGG
TGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAA
ATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAG
AAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACA
GTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCA
AAACAAAAGTAAGACCACCGCACAGCAAGCGGCCGCTGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGA
AGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCA
GAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCGT
CAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAG
GCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCT
AAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTT
GGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGC
TTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGC
AAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAG
GTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCAC
CTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCG
ATTAGTGAACGGATCGGCACTGCGTGCGCCAATTCTGCAGACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGG
GGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACA
AATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCAGTTTGGTTAATTAACCCGTGTCGGCTCC
AGATCTGGCCTCCGCGCCGGGTTTTGGCGCCTCCCGCGGGCGCCCCCCTCCTCACGGCGAGCGCTGCCACGTCAGACGAA
GGGCGCAGCGAGCGTCCTGATCCTTCCGCCCGGACGCTCAGGACAGCGGCCCGCTGCTCATAAGACTCGGCCTTAGAACC
CCAGTATCAGCAGAAGGACATTTTAGGACGGGACTTGGGTGACTCTAGGGCACTGGTTTTCTTTCCAGAGAGCGGAACAG
GCGAGGAAAAGTAGTCCCTTCTCGGCGATTCTGCGGAGGGATCTCCGTGGGGCGGTGAACGCCGATGATTATATAAGGAC
GCGCCGGGTGTGGCACAGCTAGTTCCGTCGCAGCCGGGATTTGGGTCGCGGTTCTTGTTTGTGGATCGCTGTGATCGTCA
CTTGGTGAGTAGCGGGCTGCTGGGCTGGCCGGGGCTTTCGTGGCCGCCGGGCCGCTCGGTGGGACGGAAGCGTGTGGAGA
GACCGCCAAGGGCTGTAGTCTGGGTCCGCGAGCAAGGTTGCCCTGAACTGGGGGTTGGGGGGAGCGCAGCAAAATGGCGG
CTGTTCCCGAGTCTTGAATGGAAGACGCTTGTGAGGCGGGCTGTGAGGTCGTTGAAACAAGGTGGGGGGCATGGTGGGCG
GCAAGAACCCAAGGTCTTGAGGCCTTCGCTAATGCGGGAAAGCTCTTATTCGGGTGAGATGGGCTGGGGCACCATCTGGG
GACCCTGACGTGAAGTTTGTCACTGACTGGAGAACTCGGTTTGTCGTCTGTTGCGGGGGCGGCAGTTATGGCGGTGCCGT
TGGGCAGTGCACCCGTACCTTTGGGAGCGCGCGCCCTCGTCGTGTCGTGACGTCACCCGTTCTGTTGGCTTATAATGCAG
GGTGGGGCCACCTGCCGGTAGGTGTGCGGTAGGCTTTTCTCCGTCGCAGGACGCAGGGTTCGGGCCTAGGGTAGGCTCTC
CTGAATCGACAGGCGCCGGACCTCTGGTGAGGGGAGGGATAAGTGAGGCGTCAGTTTCTTTGGTCGGTTTTATGTACCTA
TCTTCTTAAGTAGCTGAAGCTCCGGTTTTGAACTATGCGCTCGGGGTTGGCGAGTGTGTTTTGTGAAGTTTTTTAGGCAC
CTTTTGAAATGTAATCATTTGGGTCAATATGTAATTTTCAGTGTTAGACTAGTAAATTGTCCGCTAAATTCTGGCCGTTT
TTGGCTTTTTTGTTAGACGAAGCTTGGGCTGCAGGTCGACTCTAGAGGATCCCGCCACCATGGTGAGCAAGGGCGAGGAG
GTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGG
CGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACA
TCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCCGACATCCCCGATTACAAGAAGCTGTCC
TTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGTCTGGTGACCGTGACCCAGGACTCCTCCCT
GCAGGACGGCACGCTGATCTACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAGAAGAAGA
CCATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTGAAGGGCGAGATCCACCAGGCCCTGAAG
CTGAAGGACGGCGGCCACTACCTGGTGGAGTTCAAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTA
CTACGTGGACACCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCTCCGAGGGCC
GCCACCACCTGTTCCTGGGGCATGGCACCGGCAGCACCGGCAGCGGCAGCTCCGGCACCGCCTCCTCCGAGGACAACAAC
ATGGCCGTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGG
CGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCT
GGGACATCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCCGACATCCCCGATTACAAGAAG
CTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGTCTGGTGACCGTGACCCAGGACTC
CTCCCTGCAGGACGGCACGCTGATCTACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAGA
AGAAGACCATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTGAAGGGCGAGATCCACCAGGCC
CTGAAGCTGAAGGACGGCGGCCGCTACCTGGTGGAGTTCAAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGG
CTACTACTACGTGGACACCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCTCCG
AGGGCCGCCACCACCTGTTCCTGTACGGCATGGACGAGCTGTACAAGTAAGAATTCGATATCAAGCTTATCGATAATCAA
CCTCTGGATTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGC
TTTAATGCCTTTGTATCATGCTATTGCTTCCCGTATGGCTTTCATTTTCTCCTCCTTGTATAAATCCTGGTTGCTGTCTC
TTTATGAGGAGTTGTGGCCCGTTGTCAGGCAACGTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACCCCCACTGGTTGG
GGCATTGCCACCACCTGTCAGCTCCTTTCCGGGACTTTCGCTTTCCCCCTCCCTATTGCCACGGCGGAACTCATCGCCGC
CTGCCTTGCCCGCTGCTGGACAGGGGCTCGGCTGTTGGGCACTGACAATTCCGTGGTGTTGTCGGGGAAATCATCGTCCT
TTCCTTGGCTGCTCGCCTGTGTTGCCACCTGGATTCTGCGCGGGACGTCCTTCTGCTACGTCCCTTCGGCCCTCAATCCA
GCGGACCTTCCTTCCCGCGGCCTGCTGCCGGCTCTGCGGCCTCTTCCGCGTCTTCGCCTTCGCCCTCAGACGAGTCGGAT
CTCCCTTTGGGCCGCCTCCCCGCATCGATACCGTCGACCTCGAGACCTAGAAAAACATGGAGCAATCACAAGTAGCAATA
CAGCAGCTACCAATGCTGATTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTA
CCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAAT
TCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACA
CACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCAAGAGAAGGTAGAA
GAAGCCAATGAAGGAGAGAACACCCGCTTGTTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTATT
AGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGACTGTACT
EOF


# download a mouse and a human genome, and their gtfs, and combine

human_fa=GRCh38.p13.genome.fa
human_gtf=gencode.v38.basic.annotation.gtf

mouse_fa=GRCm38.p6.genome.fa
mouse_gtf=gencode.vM25.annotation.gtf

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$human_gtf".gz

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$mouse_fa".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$mouse_gtf".gz

pigz --decompress *gz

## clarify the species for each contig
sed -i "s/^>/>human_/g" "$human_fa"
sed -i "s/^>/>mouse_/g" "$mouse_fa"

grep -v "#" "$human_gtf" | sed "s/^/human_/g" > foo ; mv foo "$human_gtf"
grep -v "#" "$mouse_gtf" | sed "s/^/mouse_/g" > foo ; mv foo "$mouse_gtf"

cat $human_fa $mouse_fa captured.fa  > mixing.fa
cat $human_gtf $mouse_gtf captured.gtf  > mixing.gtf

rm $human_fa $mouse_fa $human_gtf $mouse_gtf
