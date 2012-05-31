#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;

# DEFAULTS
my $SPECIES = 'Mus musculus';
my $DIR = '/Users/remo/ANALYSIS/mmmir124';
my $result = GetOptions ('species|s=s' => \$SPECIES,
                         'dir|d=s' => \$DIR);

# ORGANIZATION
chdir($DIR);
my $OUTNAME = lc("$SPECIES\.mirbase");
$OUTNAME =~ s/ /\_/g;

download_mirbase();
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
