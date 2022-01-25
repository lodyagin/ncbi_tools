#!/usr/bin/perl
use strict;
use warnings;

# Script to convert PATRIC files to 5-column feature table files.

my $tabfile = shift or die "Must supply input filename\n";
open (my $TAB_IN, $tabfile) or die "Unable to open $tabfile\n";

my $base = $tabfile;
if ($base =~ /^(.+)\.tab$/) {
  $base = $1;
}
open (my $TBL_OUT, ">$base.tbl") or die "Unable to open feature table output file\n";

# define global variables

my $line_number = 0;
my $thisline;
my $show_first_line = 1;
my $from = 0;
my $to = 0;
my $strand = "+";

# skip first line with table headers

$thisline = <$TAB_IN>;

# main loop reads one line at a time

while ($thisline = <$TAB_IN>) {
  $thisline =~ s/\r//;
  $thisline =~ s/\n//;
  $line_number++;

  if ($thisline ne "") {
    my @items = split ('\t', $thisline);
    my $numitems = @items;
    if ($numitems > 3 && $items [3] eq "source" && $show_first_line == 1) {
      $show_first_line = 0;
      my $accn = $items [1];
      print $TBL_OUT ">Feature $accn\n";
    }

    if ($numitems > 8 && $items [3] eq "gene") {
      my $locus_tag = $items [5];
      $from = $items [6];
      $to = $items [7];
      $strand = $items [8];
      if ($strand eq "-") {
        print $TBL_OUT "$to\t$from\tgene\n";
      } else {
        print $TBL_OUT "$from\t$to\tgene\n";
      }
      print $TBL_OUT "\t\t\tlocus_tag\t$locus_tag\n";
    }

    if ($numitems > 11 && $items [3] eq "CDS") {
      my $product = $items [11];
      $from = $items [6];
      $to = $items [7];
      $strand = $items [8];
      if ($strand eq "-") {
        print $TBL_OUT "$to\t$from\tCDS\n";
      } else {
        print $TBL_OUT "$from\t$to\tCDS\n";
      }
      print $TBL_OUT "\t\t\tproduct\t$product\n";
    }
  }
}

# close input and output files

close ($TAB_IN);
close ($TBL_OUT);

