##Quick Start

###Reading Open Frame (ORF) prediction

```
perl 3UTR_orf.pl -i transcripts.fasta  -d /home/user/swissprot/swissprot -a 8 -o Test
```

###3'-UTR sequence retrieval

```
perl 3UTR_ext.pl -i1 Test_transcripts.fa -i2 Test_orfs.fasta -d /home/user/3UTR_database/3UTR.mam.fasta -a 8 -o 3UTR.fasta -x 2500 -m 20
```
