#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $usage = "\n\tUsage: perl $0 [fasta UTR] [fasta mirbase] [output]\n\n";
die $usage unless scalar(@ARGV) == 3;
die $usage unless -e $ARGV[0];
die $usage unless -e $ARGV[1];

my $UTR_FASTA = $ARGV[0];
my $MIR_FASTA = $ARGV[1];
my $OUT = $ARGV[2];

open(OUT,">$OUT");
print OUT "transcript_id\tmiRNA_id\t7mer\-A1\t7mer-m8\t8mer\tTOT\n";

my $UTR = Bio::SeqIO->new(-file => $UTR_FASTA,
                          -format => 'fasta');

my $MIR = Bio::SeqIO->new(-file => $MIR_FASTA,
                          -format => 'fasta');

my $mirna_seqs = get_mirna();

while(my $seq = $UTR->next_seq) {
  my $tid = $seq->id;
  foreach my $id(keys %$mirna_seqs) {
    my $counts = scan_utr($seq->seq,$mirna_seqs->{$id});
    my $sum = $counts->{'7mer-A1'}->{count}+$counts->{'7mer-m8'}->{count}+$counts->{'8mer'}->{count};
    print OUT join("\t",($tid,$id,$counts->{'7mer-A1'}->{count},$counts->{'7mer-m8'}->{count},$counts->{'8mer'}->{count},$sum));
    print OUT "\n";
  }
}

close(OUT);

exit;

sub get_mirna {
  my $href = {};
  while(my $seq = $MIR->next_seq) {
    $href->{$seq->id} = $seq->seq;
  }
  return $href;
}

sub scan_utr {

	my $utr = shift;
	my $mirna = shift;
	my $motif = '';
	my $motifs = {};

	# CHANGE 'U' WITH 'T'
	$mirna =~ tr/uU/tT/;

	# EXTRACT NUCLEOTIDES 2..8 FROM MIRNA
	my @bases = split(//,$mirna);
	my $core = join('',@bases[1..7]);
	my $seed = join('',@bases[1..6]);

	my $core_rev = reverse($core);
	$core_rev =~ tr/ACGTacgt/TGCAtgca/;
	my $seed_rev = reverse($seed);
	$seed_rev =~ tr/ACGTacgt/TGCAtgca/;
	
	# 7mer-A1
	# NUCLEOTIDES 2..7 INVERTED + 'A' AT THE END
	$motifs->{'7mer-A1'}->{count} = 0;
	$motif = "$seed_rev".'A';
	$motifs->{'7mer-A1'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'7mer-A1'}->{count} ++;
	}

	# 7mer-m8
	# NUCLEOTIDES 2..8 INVERTED
	$motifs->{'7mer-m8'}->{count} = 0;
	$motif = "$core_rev";
	$motifs->{'7mer-m8'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'7mer-m8'}->{count} ++;
	}

	# 8mer
	# NUCLEOTIDES 2..8 INVERTED + 'A' AT THE END
	$motifs->{'8mer'}->{count} = 0;
	$motif = "$core_rev".'A';
	$motifs->{'8mer'}->{motif} = $motif;
	while($utr =~ /$motif/ig) {
		$motifs->{'8mer'}->{count} ++;
	}

	for(my $c=1; $c<=$motifs->{'8mer'}->{count}; $c++) {
		$motifs->{'7mer-A1'}->{count} --;
		$motifs->{'7mer-m8'}->{count} --;
	}

	return $motifs;
}
