##Installation

- created 15/02/2017

##Requirements


 - Perl http://www.perl.org/
 - Bioperl http://www.bioperl.org/
 - ncbi-blast ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
 - Swissprot ftp://ftp.ncbi.nlm.nih.gov/blast/db/
 - 3UTR.mam.fasta http://utrdb.ba.itb.cnr.it/home/download


##Installation

The instructions for the installation are based on Ubuntu Linux system.

**1. Install the dependent tools and databases.**

- Perl 5.6.1 or higher

- Bioperl

```
[user@ubuntu:~]$ git clone https://github.com/bioperl/bioperl-live.git
[user@ubuntu:~]$ cd bioperl-live
[user@ubuntu:~]$ perl Build.PL
[user@ubuntu:~]$ ./Build test
[user@ubuntu:~]$ ./Build install
```

- ncbi-blast

```
[user@ubuntu:~]$ sudo apt-get install ncbi-blast+
```

- Swissprot

The Swissprot database is pre-formatted which is ready to use with BLAST.


- 3UTR.mam.fasta 

This is optional, and you can use other 3'-UTR databases depending on your species.

Due to the duplicated sequence ids in 3UTR.mam.fasta database, errors will be reported during formatting the database. Please either remove the sequences with duplicated ids or use the reformatted version of 3UTR.mam.fasta database which can be downloaded in this project in github.


**2. Install the ExUTR.**

```
[user@ubuntu:~]$ chmod 777 ExUTR/*
[user@ubuntu:~]$ perl 3UTR_orf.pl

3UTR_orf.pl v1.0 15-02-2017


Usage:


perl 3UTR_orf.pl [ -i input_file_path ] [ -d swissprot_database_path ] [ -a blast_thread ] [ -o output_name ]


Description:

    Parameters:
    
    -i input file (transcripts) path
    -d directory path to Swissprot database
    -a the number of threads for blast search, default: 2.
    -o prefix name for output files, default: project

Example: perl 3UTR_orf.pl -i transcripts.fa -d path_to_swissprot -a 8 -o project

```
```
[user@ubuntu:~]$ perl 3UTR_ext.pl

3UTR_ext.pl v1.0 15-02-2017


Usage:


perl 3UTR_ext.pl [ -i1 input_transcripts ] [ -i2 input_orfs ] [ -a blast_thread ] [ -o output_name.fa ] [ -x maximum 3UTR length ] [ -m minimum 3UTR length ] [ -d 3UTR_database ]


Description:

    Parameters:
    
    -i1 input transcripts FASTA file
    -i2 input ORF FASTA file
    -a the number of threads for blast search, default: 2.
    -o prefix name for output files, default: 3UTR.fasta
    -x maximum 3UTR length allowed, default: 2000.
    -m minimum 3UTR length allowed, default: 20.
    -d directory path to the 3UTR database


Example: perl 3UTR_ext.pl -i1 transcripts.fa -i2 orf.fa -d path_to_3UTR_database -a 8 -o 3UTR.fasta -x 2500 -m 30

Email bugs to zixia.huang@ucdconnect.ie

```
