#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# DEFINITION OF THE FIELDS AND THEIR ORDER IN THE FASTA
# DOWNLOADED FROM BIOMART. TO BE EDITED ACCORDINGLY.
my @info = qw( geneid
               transcriptid
               symbol
               seqregion
               biotype
               strand
               start
               end );

# DEFINITION OF THE FIELD UNIQUE IN THE FASTA.
# TO BE EDITED ACCORDINGLY.
my $unique = 'transcriptid';

# USAGE AND CHECK OF THE PARAMETERS.
my $usage = "\n\tUsage: perl $0 [fasta from biomart] [length cutoff] [coding filter (1 or 0)]\n\n";
die $usage unless scalar(@ARGV) == 3;
die $usage unless -e $ARGV[0];
die $usage unless $ARGV[1] >= 1;

my $file = $ARGV[0];
my $lcutoff = $ARGV[1];
my $coding = $ARGV[2];

my $seqio = Bio::SeqIO->new(-file => $file,
                            -format => 'fasta');

my $out = Bio::SeqIO->new(-file => ">mod_$file",
                          -format => 'fasta');

my $anno = "anno_$file";
$anno =~ s/\.fasta//;
$anno =~ s/\.fa//;
$anno =~ s/\.txt//;
$anno .= '.txt';
open(OUT,">$anno");
print OUT $unique;
for(my $c = 0; $c <= $#info; $c++) {
  my $info = $info[$c];
  next if $info eq $unique;
  print OUT "\t$info";
}
print OUT "\n";

while(my $seq = $seqio->next_seq) {
  next unless $seq->length >= $lcutoff;
  next if($coding && $seq->id !~ /protein_coding/);
  my $lid = $seq->id;
  my @val = split(/\|/,$lid);
  my $string;
  for(my $c = 0; $c <= $#info; $c++) {
    my $info = $info[$c];
    if($info eq $unique) {
      $seq->id($val[$c]);
    }
    else {
      my $val = $val[$c] || 'NA';
      $string .= "\t$val";
    }
  }
  print OUT $seq->id.$string."\n";
  $out->write_seq($seq);
}
