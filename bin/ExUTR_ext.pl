#!/usr/bin/env perl 


=head1 NAME
3UTR_ext.pl
=head1 DESCRIPTION
This script allows you to extract 3-UTR region from the assembled transcripts.
     
 Author : Zixia Huang
 Email : zixia.huang@ucdconnect.ie
 Date : 16-08-2017
 Version : 0.1.0
=cut


use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;


my $usage = q/
3UTR_ext.pl v1.0 16-08-2017
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
/;


die $usage."\n" if (@ARGV < 1);


## Set up input, output files and other arguments.


my $input_tr;
my $input_orf;
my $length_max = 2000;
my $length_min = 20;
my $thread = 2;
my $output = '3UTR.fasta';
my $database;




GetOptions (

      'i1:s' => \$input_tr,
      'i2:s' => \$input_orf,
      'a:n' => \$thread,
      'o:s' => \$output,
      'x:n' => \$length_max,
      'm:n' => \$length_min,
      'd:s' => \$database,

       );



###############################################################
#          BLAST transcripts to corresponding ORFs            #
###############################################################



print "\n";
print "********** WELCOME **********\n\n";


## Check the input files.
open(TA,$input_tr) or die ("$!, stopped");
open(TB,$input_orf) or die ("$!, stopped");


print "BLAST-ing transcripts to predicted ORFs......\n";


`makeblastdb -in $input_orf -hash_index -parse_seqids -dbtype prot`;
`blastx -query $input_tr -db $input_orf -out tmp_blast_report -num_alignments 1 -num_descriptions 1 -num_threads $thread -evalue 0.000001 2>error.log`;


###############################################################
#           Calculate 3UTR position in transcripts            #
###############################################################


print "\nFinish BLAST search and start calculating 3UTR positions......\n";



open(AB,'>tmp_UTR_forward');
open(CD,'>tmp_UTR_reverse');
open(EF,'>tmp_frame');



my $crit = 0.000001;
my $br = new Bio::SearchIO (-file => 'tmp_blast_report', format => 'blast');

                RESULTS: while (my $result = $br -> next_result)
                        
                        {
                            my $query_name = $result -> query_accession;

                            my $query_length = $result -> query_length;
 
                            my $query_desc = $result -> query_description;
                       
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
 
                                                    if($query_name eq $hid) 
  
                                                    {
                                               
                                                    my $frame = ($hsp -> query -> frame+1)*($hsp -> query -> strand);
                                                    my $start = $hsp -> start('query');
                                                    my $end = $hsp -> end('query');


                                                    print EF $frame."\n";

                                                    # Process the transcripts with 5'-oriented
                                                    if ($frame > 0 and ($query_length - $end) >= $length_min 
                                                        and ($query_length - $end) <= $length_max) 

                                                      { 

                                                          print AB $query_name."\t".$end."\n";

                                                      }
                   
                                                     # Process the transcripts with 3'-oriented
                                                     if ($frame < 0 and ($query_length - $end) >= $length_min 
                                                        and ($query_length - $end) <= $length_max) 

                                                      { 

                                                          print CD $query_name."\t".$end."\n";

                                                      }
                                                 
                                                  }
                                              
                                              }
                                             
                                          }
                                           
                                      }
                                  
                                 }
                    


               
###############################################################
#             Extract 3UTR region from transcripts            #
###############################################################



print "\nStart to extract 3UTR sequences from transcripts......\n";
`ls`;


## This subroutine is used to extract 3-UTR regions from transcripts based on the positions of the stop codons.


sub EXTRACT_3UTR {
  
   my ($a,$b,$c) = @_;

   my %NR;

   open($a,"$a");
   open($c,">$c");

   my @file = <$a>;

   foreach my $t (@file) 

   { 
      chomp($t);

      my @line = split /\t/, $t;

      $NR{$line[0]} = $line[1];

   }
  
   my $in = Bio::SeqIO -> new (-file => $input_tr);

       while (my $seqio = $in -> next_seq)

       {
           my $name = $seqio -> primary_id;

           my $length = $seqio -> length;

           my $desc = $seqio -> desc;

           if (defined($NR{$name}))

              {

                 if ($b eq 'A')  # Extract 3UTR from 5-oriented transcripts;

                 {

                   my $UTR = $seqio -> subseq($NR{$name}, $length);

                   print $c '>'.$name."\t".$desc."\n".$UTR."\n";

                 }

                 if ($b eq 'B')  # Extract 3UTR from 3-oriented transcripts;

                 {

                   my $UTR = $seqio -> revcom -> subseq($NR{$name}, $length);

                   print $c '>'.$name."\t".$desc."\n".$UTR."\n"; 
                 
                 }


              }

          }

     }



my @forward = qw /tmp_UTR_forward A tmp_UTR_forward.fa/;        

EXTRACT_3UTR (@forward);
                 
        
my @reverse = qw /tmp_UTR_reverse B tmp_UTR_reverse.fa/;

EXTRACT_3UTR (@reverse);



###############################################################
#                    Output 3UTR sequences                    #
###############################################################


`cat tmp_UTR_forward.fa tmp_UTR_reverse.fa >$output`;


print "\nGenerating 3UTR file: $output.\n";



###############################################################
#                       3UTR validation                       #
###############################################################


`makeblastdb -in $database -hash_index -parse_seqids -dbtype nucl 2>error.log`;
`blastn -query $output -db $database -out tmp_blast_validation -num_alignments 1 -num_descriptions 1 -num_threads $thread -evalue 1e-6 2>error.log`;



 open(AA,'>orf_validation'); ## Summary of predicted 3'-UTR which had BLAST hits in UTRdb.


 print AA 'Predicted 3UTR'."\t".'Hit in 3UTR Database.'."\t".'3UTR length'."\t".'Gene name'."\n";

 my $crit_v = 0.000001;
 my $count=0;

 my $val = new Bio::SearchIO -> new (-file => 'tmp_blast_validation', format => 'blast');

                RESULTS: while (my $result_v = $val -> next_result)

                  {
                       my $query_name_v = $result_v -> query_accession;
                       my $query_length_v = $result_v -> query_length;
                       my $query_desc_v = $result_v -> query_description;

                          HITS: while (my $hit_v = $result_v -> next_hit)

                            {
                                my $hobs_v = $hit_v -> significance;

                                if ($hobs_v > $crit_v) 

                                     {
                                       last HITS;
                                                  }

                                HSPS: while (my $hsp_v = $hit_v -> next_hsp)
                                      
                                    {
                                          my $eobs_v = $hsp_v -> evalue;

                                          if ($eobs_v > $crit_v) 
                                     
                                             {
                                                last HSPS;
                                                            }

                                          if ($eobs_v <= $crit_v)

                                             {
                                                my $hid_v = $hit_v -> accession;

                                                $count++;

                                                print AA $query_name_v."\t".$hid_v."\t".$query_length_v."\t".$query_desc_v."\n";

                                                last;         
                                                            }

                                        last;      
                                                   }

                                 last;      
                                          }


                                 }

  `rm tmp*`;
  `rm *phr *pin *psq *phd *phi *pog *psd *psi`;


  print "\n********SUMMARY********\n";
  print "\nThe number of transcripts with predicted 3UTR : \n";

  print `grep -c '>' $output`;

  print "\n";

  print "The number of predicted 3UTR with hits in UTRdb database : \n";
  print "$count\n";
  
