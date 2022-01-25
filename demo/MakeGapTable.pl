#!/usr/bin/perl

use strict;

my $infile;
my $outfile;
my $thisline;
my $id;
my $new_id;
my $sequence;
my $mingap;

$infile = shift || die "Must specify input filename!\n";
$outfile = shift;
if (!$outfile) {
  $outfile = $infile;
  $outfile =~ s/\.fsa$//;
  $outfile .= ".tbl";
}

$mingap = shift;
if ($mingap < 1) {
  $mingap = 1;
}
 
open (IN, $infile) || die "Unable to open $infile!\n";
open (OUT, ">$outfile") || die "Unable to open $outfile!\n";
while ($thisline = <IN>) {
  $thisline =~ s/\r//g;
  $thisline =~ s/\n//g;
  if ($thisline =~ /\>\s*(\S+)/) {
    $new_id = $1;
    if ($id) {
      &process_seq($id, $sequence, $mingap);
    }
    $id = $new_id;
    $sequence = "";
  } else {
    $sequence .= $thisline;
  }
}
if ($id) {
  &process_seq($id, $sequence);
} else {
  print "No sequences found!\n";
}
close(IN);
close(OUT);

sub process_seq
{
  my $id = shift(@_);
  my $sequence = shift(@_);
  my $mingap = shift(@_);
  my $pos;
  my $offset = 0;
  my $len;

  print OUT ">Feature $id\n";
  $sequence =~ s/\s//g;
  $sequence = uc($sequence);
  $pos = index($sequence, "N");
  while ($pos != -1) {
    $offset += $pos;
    $sequence = substr($sequence, $pos);
    $len = 0;
    while(substr($sequence, $len, 1) eq "N") {;
      $len++;
    }
    if ($len >= $mingap) {
      print OUT $offset + 1;
      print OUT "\t";
      print OUT $offset + $len;
      print OUT "\tassembly_gap\n";
      print OUT "\t\t\tgap_type\twithin scaffold\n";
      print OUT "\t\t\tlinkage_evidence\tpaired-ends\n";
      print OUT "\t\t\testimated_length\t";
      print OUT $len;
      print OUT "\n";
    }
    $sequence = substr($sequence, $len);
    $offset += $len;
    $pos = index($sequence, "N");
  }
}

