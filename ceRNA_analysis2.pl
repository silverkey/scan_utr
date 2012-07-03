#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# TO BE RUN IN A UNIQUE FOLDER CONTAINING ALL THE INPUT FILES!
my $usage = "\n\tTO BE RUN IN A UNIQUE FOLDER CONTAINING ALL THE INPUT FILES!\n";
$usage .= "\n\tUsage: perl $0 [fasta NC] [fasta UTR] [fasta mirbase] [UTR annotation or 'NA'] [field from anno or 'NA'] [min NC sites] [min UTR sites]\n\n";
die $usage unless scalar(@ARGV) == 7;
die "Cannot find file $ARGV[0] \n\n $usage" unless -e $ARGV[0];
die "Cannot find file $ARGV[1] \n\n $usage" unless -e $ARGV[1];
die "Cannot find file $ARGV[2] \n\n $usage" unless -e $ARGV[2];
if($ARGV[3] ne 'NA'){die "Cannot find file $ARGV[3] \n\n $usage" unless -e $ARGV[3]}
die "Wrong parameter $ARGV[5] \n\n $usage" unless $ARGV[5] =~ /^\d+$/;
die "Wrong parameter $ARGV[6] \n\n $usage" unless $ARGV[6] =~ /^\d+$/;

my $NC_FASTA = $ARGV[0];
my $UTR_FASTA = $ARGV[1];
my $MIR_FASTA = $ARGV[2];
my $UTR_ANNO = $ARGV[3];
my $ANNO_FIELD = $ARGV[4];
my $MIN_NC_SITE = $ARGV[5];
my $MIN_UTR_SITE = $ARGV[6];

my $MIRNA_SEQS = get_mirna();

my $ANNO = get_anno() if $UTR_ANNO ne 'NA';

my $NC_SITES = site_search($NC_FASTA,$MIN_NC_SITE);
my $UTR_SITES = site_search($UTR_FASTA,$MIN_UTR_SITE);

open(SUM,">ceRNA_summary.xls");
print SUM "noncoding\tmiRNA\tsites\n";
foreach my $tid(keys %{$NC_SITES->{tid}}) {
  my $nc = $NC_SITES->{tid}->{$tid};
  foreach my $mir(keys %$nc) {
    my $n = $nc->{$mir};
    print SUM join("\t",$tid,$mir,$n)."\n";
    my @target = keys(%{$UTR_SITES->{mid}->{$mir}});
    next unless scalar(@target) >= 1;
    open(OUT,">$tid\__$mir\__$n.xls");
    print OUT "transcript\tannotation\tsites\n";
    foreach my $target(@target) {
      my $anno = $ANNO->{$target} || 'NA';
      print OUT join("\t",$target,$anno,$UTR_SITES->{mid}->{$mir}->{$target})."\n";
    }
    close(OUT);
  }
}

sub site_search {
  my $fasta = shift;
  my $min = shift;
  my $href = {};
  my $target = Bio::SeqIO->new(-file => $fasta,
                               -format => 'fasta');
  my $out = "$fasta\_$MIR_FASTA";
  open(OUT,">$out");
  print OUT "transcript_id\tmiRNA_id\t7mer\-A1\t7mer-m8\t8mer\tTOT\n";
  while(my $seq = $target->next_seq) {
    my $tid = $seq->id;
    foreach my $id(keys %$MIRNA_SEQS) {
      my $counts = scan_utr($seq->seq,$MIRNA_SEQS->{$id});
      my $sum = $counts->{'7mer-A1'}->{count}+$counts->{'7mer-m8'}->{count}+$counts->{'8mer'}->{count};
      next unless $sum >= 1;
      print OUT join("\t",($tid,$id,$counts->{'7mer-A1'}->{count},$counts->{'7mer-m8'}->{count},$counts->{'8mer'}->{count},$sum));
      print OUT "\n";
      if($sum >= $min) {
        $href->{tid}->{$tid}->{$id} = $sum;
        $href->{mid}->{$id}->{$tid} = $sum;
      }
    }
  }
  close(OUT);
  return $href;
}

sub get_anno {
  my $href = {};
  open(IN,$UTR_ANNO);
  my $head = <IN>;
  chomp($head);
  my @head = split(/\t/,$head);
  while(my $row = <IN>) {
    chomp($row);
    my @field = split(/\t/,$row);
    for(my $c = 0; $c <= $#head; $c++) {
      next unless $head[$c] eq $ANNO_FIELD;
      $href->{$field[0]} = $field[$c];
    }
  }
  return $href;
}

sub get_mirna {
  my $mir = Bio::SeqIO->new(-file => $MIR_FASTA,
                            -format => 'fasta');
  my $href = {};
  while(my $seq = $mir->next_seq) {
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
__END__

