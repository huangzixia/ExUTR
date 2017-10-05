#!/usr/bin/perl


=head1 NAME
3UTR_orf.pl
=head1 DESCRIPTION
This script allows you to predict potential ORFs from the assembled transcripts.
     
 Author : Zixia Huang
 Email : zixia.huang@ucdconnect.ie
 Date : 16-08-2017
 Version : 0.1.0
=cut


use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

my $usage = q/
3UTR_orf.pl v1.0 15-02-2017
Usage:
perl 3UTR_orf.pl [ -i input_file_path ] [ -d swissprot_database_path ] [ -a blast_thread ] [ -o output_name ] [ -l lib_type ]
Description:
    Parameters:
    
    -i input file (transcripts) path
    -d directory path to Swissprot database
    -a the number of threads for blast search, default: 2.
    -o prefix name for output files, default: project
    -l library type (un: unstranded; ff: first-forward; fr: first-reverse)
Example: perl 3UTR_orf.pl -i transcripts.fa -d path_to_swissprot -a 8 -o project -l un
/;


die $usage."\n" if (@ARGV < 1);



##Setting input, output files and other arguments.


my $file_input;
my $crit = 0.000001;
my $database = "swissprot"; 
my $thread = 2;
my $output = "project";
my $lib_type;


GetOptions (

      'i:s' => \$file_input,
      'd:s' => \$database,
      'a:n' => \$thread,
      'o:s' => \$output,
      'l:s' => \$lib_type,

       );


open(TA,$file_input) or die ("$!, stopped");
open(TB,$database.'.pal') or die ("$!, stopped");



###############################################################
#                  Potential ORF Translation                  #
###############################################################


print "\n";
print "********** WELCOME **********\n";
print "\n";
print "Start time:\n";
print `date`;
print "\n";
print "Starting to translate potential ORFs......\n";




 ##copy and reformat the input FASTA file
`sed 's/>/>lcl|/' $file_input >tmp_rearrange.fa`; 



my $file_input_1='tmp_rearrange.fa';
my @length;
my @sort_length;
my %NR;


open(AA,'>tmp_all_transcripts_id');
open(OA,'>tmp_orf1.fa');
open(OB,'>tmp_orf2.fa');


my $in = Bio::SeqIO -> new ( -file => $file_input_1 ); ## Read input FASTA file

     while (my $seqio = $in -> next_seq) 
        
       {

           my $name = $seqio -> primary_id;
           my $seq = $seqio -> seq;
           my $seq_length = $seqio -> length;

           if ($lib_type eq 'un') 

           {

           print AA $name,"\n";   ## Write all the transcript ids into 'AA'.

           ## Each transcript is tranlated in 6 frames with both forward and reverse directions.
           ## For each direction the longest ORF will be extracted.

           my $orf1 = $seqio -> translate (-orf => 'longest') -> seq;

  
           my $orf_re1 = $seqio -> revcom -> translate (-orf => 'longest') -> seq;

           

           $length[0] = length($orf1);

           $length[1] = length($orf_re1);

  
           $NR{$length[0]} = $orf1;

           $NR{$length[1]} = $orf_re1;


           ## For each transcript, sort the potential ORFs in length.
           ## Then save them in respective files. 

           @sort_length = sort { $b <=> $a } @length;
 
           print OA '>'.$name."\n".$NR{$sort_length[0]}."\n";  ## Write the longer ORF into 'OA';
           print OB '>'.$name."\n".$NR{$sort_length[1]}."\n";  ## Write the shorter ORF into 'OB';

           }


           if ($lib_type eq 'ff')

           {

           print AA $name,"\n";   ## Write all the transcript ids into 'AA'.

           my $orf1 = $seqio -> translate (-orf => 'longest') -> seq;

           print OA '>'.$name."\n".$orf1."\n";

           }


           if ($lib_type eq 'fr')

           {

           print AA $name,"\n";   ## Write all the transcript ids into 'AA'.

           my $orf_re1 = $seqio -> revcom -> translate (-orf => 'longest') -> seq;

           print OA '>'.$name."\n".$orf_re1."\n";

           }
      
       }


print "\n";
print "Finish ORF translation.\n";
print "\n-----\nDone\n-----\n\n";
print "Now Blasting the translated ORFs to Swissprot......\n";



###############################################################
#                  Analysing blast reports                    #
###############################################################



## Blast the longer ORFs to Swissprot.


`blastp -query tmp_orf1.fa -db $database -out tmp_blast_orf1 -num_descriptions 1 -num_alignments 1 -num_threads $thread -evalue 0.000001 2>error.log`;


## Analyzing the blast report 'tmp_blast_orf1'.
    

my @file1 = qw /tmp_blast_orf1 tmp_blast_hit_1/;


## The subroutine is used to parse blast report, looking for ORFs which have BLAST hits in Swissprot.


sub PARSE_BLAST_REPORT {

         my ($a,$b) = @_;

         open ($b,">$b");

         my $br = new Bio::SearchIO -> new (-file => @_, format => 'blast');

                RESULTS: while (my $result = $br -> next_result)

                  {
                       my $query_name = $result -> query_accession;

                          HITS: while (my $hit = $result -> next_hit)

                            {
                                my $hobs = $hit -> significance;

                                if ($hobs > $crit) 

                                     {
                                       last HITS;
                                                  }

                                HSPS: while (my $hsp = $hit -> next_hsp)
                                      
                                    {
                                          my $eobs = $hsp -> evalue;

                                          if ($eobs > $crit) 
                                     
                                             {
                                                last HSPS;
                                                            }

                                          if ($eobs <= $crit)

                                             {
                                                my $hid = $hit -> accession;
                                                
                                                print $b $query_name."\t".$hid."\n";

                                                last;         
                                                            }

                                        last;      
                                                   }

                                 last;      
                                          }


                                 }


                        }



PARSE_BLAST_REPORT (@file1);

## Extract those ORFs without BLAST hits.

            `cut -f 1 tmp_blast_hit_1 | sort >tmp_blast_hit_1_sort`;

            `sort tmp_all_transcripts_id >tmp_all_transcripts_id_sort`;

            `comm -3 tmp_all_transcripts_id_sort tmp_blast_hit_1_sort >tmp_no_hit_1`;

            `makeblastdb -in tmp_orf2.fa -hash_index -parse_seqids -dbtype prot 2>error.log`;

            `blastdbcmd -entry_batch tmp_no_hit_1 -db tmp_orf2.fa -dbtype prot -out tmp_no_hit_orf1.fa 2>error.log`;



## For those without blast hits, blast their corresponding shorter ORFs to Swissprot.


`blastp -query tmp_no_hit_orf1.fa -db $database -out tmp_blast_orf2 -num_alignments 1 -num_descriptions 1 -num_threads $thread -evalue 0.000001 2>>error.log`;


my @file2 = qw /tmp_blast_orf2 tmp_blast_hit_2/;


## Analize the shorter ORFs with BLAST hits in Swissprot.


PARSE_BLAST_REPORT (@file2);


`cut -f 1 tmp_blast_hit_2 | sort >tmp_blast_hit_2_sort`;



## Finish blast report interpretation.

print "Blast and report analysis.\n";
print "\n-----\nDone\n-----\n\n";



###############################################################
#                        Final Output                         #
###############################################################


## Generate a FASTA file of transcripts with potential ORFs in Swissprot.
## Generate a FASTA file of corresponding ORFs of those transcripts.


            `cat tmp_blast_hit_1 tmp_blast_hit_2 >tmp_blast_hit`;

            `cat tmp_blast_hit_1_sort tmp_blast_hit_2_sort >tmp_blast_hit_name`;

            `makeblastdb -in $file_input_1 -hash_index -parse_seqids -dbtype prot 2>error.log`;

            `cut -f 2 -d '|' tmp_blast_hit_name >tmp_blast_hit_name_format`;

            `blastdbcmd -entry_batch tmp_blast_hit_name_format -db $file_input_1 -dbtype prot -out tmp_transcripts_orf.fa 2>error.log`;



## Annotate these transcripts (Swissprot ID will be shown on the sequence title)


print "Annotating the transcripts......\n\n";

my %NO;

my $to = Bio::SeqIO -> new (-file => 'tmp_transcripts_orf.fa');
  
         while (my $toi=$to->next_seq) 
       
          {
                my $toname = $toi -> primary_id;
                my $toseq = $toi -> seq;
                
                $NO{$toname}=$toseq;
           }


open(CC,'tmp_blast_hit');

open(OP, ">tmp_transcripts_with_orf.fa");

          my @lines = <CC>;
    
          foreach my $t (@lines) 

           {
             chomp($t);
             my @a = split /\t/,$t;
    
                   if (defined $NO{$a[0]})  

                      {

                          print OP '>'.$a[0]." ".$a[1]."\n".$NO{$a[0]}."\n";

                      }
 
           }


## Calculate the size of these files containing the transcript ids which have BLAST hits. If empty, it will be excluded from further analysis.


   my $file_size_1= -s 'tmp_blast_hit_1_sort';
   my $file_size_2= -s 'tmp_blast_hit_2_sort';


`makeblastdb -in tmp_orf1.fa -hash_index -parse_seqids -dbtype prot 2>error.log`;


print "Extracting ORF sequences with the stop codon......\n\n";


## This subroutine allows to extract the ORFs with stop codon.



sub EXTRACT_ORF {

    my ($c,$d,$e) = @_;

    my %RR;

    open($c,"$c");
    open($e,">$e");
    
    my $in_orf = Bio::SeqIO -> new (-file => $d);

        while (my $seq_orf = $in_orf -> next_seq) 

            { 
               my $name_orf = $seq_orf -> primary_id;
               my $sequence_orf = $seq_orf -> seq;

               if ($sequence_orf =~ /\*/) 

                   {
                       $RR{$name_orf} = $sequence_orf;
                   }

             }

 

    my @orf_lines = <$c>;

    foreach my $line (@orf_lines)

             {  
                chomp($line);

                if (defined $RR{$line}) 

                   {
                      print $e '>'.$line."\n".$RR{$line}."\n";
                   }

              }

         }



if ($file_size_1 > 0 )

   {  
      my @file_orf_1= qw / tmp_blast_hit_1_sort tmp_orf1.fa tmp_orf1_sp.fa /;
      EXTRACT_ORF (@file_orf_1);
   }

if ($file_size_2 > 0 )

   {
      my @file_orf_2= qw / tmp_blast_hit_2_sort tmp_orf2.fa tmp_orf2_sp.fa /;
      EXTRACT_ORF (@file_orf_2);
   }




`cat *sp.fa >"$output"_orfs.fa`; ## merge the predicted orf sequences

`makeblastdb -in tmp_transcripts_with_orf.fa -hash_index -parse_seqids -dbtype nucl 2>error.log`;



open(OQ,"$output''_orfs.fa");

open(OZ,'>tmp_transcripts_with_sc'); ## sc:stop codon



my $n = $output.'_orfs.fa';

my $sc = Bio::SeqIO -> new (-file => $n);

  while (my $seq_sc = $sc -> next_seq )

    {
       my $name_sc = $seq_sc -> primary_id;

       print OZ $name_sc."\n";

    }


## Extract the corresponding transcripts of the predicted ORFs.

`blastdbcmd -entry_batch tmp_transcripts_with_sc -db tmp_transcripts_with_orf.fa -dbtype nucl -out "$output"_transcripts.fa 2>error.log`;

     
## Remove the intermediate files.
`rm tmp*`; 


print "Finish time:\n";
print `date`;
print "\n";
print "Excellent!! Everything Done !! Now run 3UTR_ext.pl";
print "\n";

