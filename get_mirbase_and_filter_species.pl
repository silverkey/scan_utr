#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\n\tUsage: perl $0 [species] [download mirbase? (1 or 0)]\n\n";
die $usage unless scalar(@ARGV) == 2;
die $usage unless $ARGV[0] =~ /^[A-Z][a-z]+ [a-z]+$/;

my $SPECIES = $ARGV[0];
my $DOWNLOAD = $ARGV[2];

my $OUTNAME = lc("$SPECIES\.mirbase");
$OUTNAME =~ s/ /\_/g;

download_mirbase() if $DOWNLOAD;
parse_mirbase();

sub download_mirbase {
  my $mirbase = 'ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz';
  my $cmd = "wget $mirbase";
  system($cmd);
  system('gunzip mature.fa.gz');
}

sub parse_mirbase {
  my $mirbase = Bio::SeqIO->new(-file => 'mature.fa',
                                -format => 'fasta');

  my $filtered = Bio::SeqIO->new(-file => ">$OUTNAME",
                                 -format => 'fasta');

  while(my $seq = $mirbase->next_seq) {
    $filtered->write_seq($seq) if $seq->desc =~ /$SPECIES/;
  }
}
