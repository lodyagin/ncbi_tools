#!/usr/bin/perl

# ftp ftp://ftp.expasy.org/databases/enzyme/enzyme.dat
# ftp ftp://ftp.expasy.org/databases/enzyme/enzclass.txt

# ./get_ec.pl enzyme.dat enzclass.txt

# /am/ncbiapdata/scripts/misc/txt2inc.sh ecnum_ambiguous.txt 
# /am/ncbiapdata/scripts/misc/txt2inc.sh ecnum_deleted.txt 
# /am/ncbiapdata/scripts/misc/txt2inc.sh ecnum_replaced.txt 
# /am/ncbiapdata/scripts/misc/txt2inc.sh ecnum_specific.txt 

$enzfile = shift || die "Must supply enzyme.dat filename\n";
$clsfile = shift || die "Must supply enzclass.txt filename\n";
open (ENZIN, $enzfile) || die "Unable to open $enzfile\n";
open (CLSIN, $clsfile) || die "Unable to open $clsfile\n";
open (OUT_LIVE, ">ecnum_specific.txt") || die "Unable to open ecnum_specific.txt\n";
open (OUT_DEL, ">ecnum_deleted.txt") || die "Unable to open ecnum_deleted.txt\n";
open (OUT_TRANS, ">ecnum_replaced.txt") || die "Unable to open ecnum_replaced.txt\n";
open (OUT_AMB, ">ecnum_ambiguous.txt") || die "Unable to open ecnum_ambiguous.txt\n";

while ($thisline = <ENZIN>) {
  $thisline =~ s/\r//;
  $thisline =~ s/\n//;
  if ($thisline =~ /^CC/) {
    #ignore comment lines
  } elsif ($thisline =~ /^ID\s+(.*)/) {
    $current_id = $1;
    #by default, entry type is 1
    $entry_type = 1;
    $print_id = 1;
    $add_space = 0;
  } elsif ($thisline =~ /^DE\s+(.*)/) {
    $disposition = $1;
    if ($disposition =~ /Deleted entry/) {
      print OUT_DEL "$current_id";
      $entry_type = 2;
    } elsif ($disposition =~ /Transferred entry: (.*)/) {
      print OUT_TRANS "$current_id";
      $entry_type = 3;
      $disposition = $1;
    }
    if ($entry_type == 1) {
      if ($print_id == 1) {
        print OUT_LIVE "$current_id\t";
        $print_id = 0;
      }
      if ($add_space == 1) {
        print OUT_LIVE " ";
      }
      #use substitution to remove trailing period
      $disposition =~ s/\.\s*$//;
      print OUT_LIVE "$disposition";
      if ($disposition !~ /-\s*$/) {
        $add_space = 1;
      }
    } elsif ($entry_type == 3) {
      $next_id = $disposition;
      #use substitution to remove and
      $next_id =~ s/ and//;
      $next_id =~ s/and //;
      #use substitution to remove commas (note use of g for global)
      $next_id =~ s/,//g;
      #use substitution to remove trailing period
      $next_id =~ s/\.\s*$//;
      #use substitution to replace spaces with tabs (note use of g for global)
      $next_id =~ s/\s+/\t/g;
      print OUT_TRANS "\t$next_id";
    }
  } elsif ($thisline =~ /^\/\//) {
    if ($entry_type == 1) {
      print OUT_LIVE "\n";
    } elsif ($entry_type == 2) {
      print OUT_DEL "\n";
    } elsif ($entry_type == 3) {
      print OUT_TRANS "\n";
    }
  }
}

while ($thisline = <CLSIN>) {
  $thisline =~ s/\r//;
  $thisline =~ s/\n//;
  $thisline =~ s/\.\s+/\./;
  if ($thisline =~ /^(.+- )\s+(.*)/) {
    $ec_num = $1;
    $ec_name = $2;
    #use substitution to delete spaces (note use of g for global)
    $ec_num =~ s/\s+//g;
    #use substitution to remove trailing period
    $ec_name =~ s/\.\s*$//;
    print OUT_AMB "$ec_num\t$ec_name\n";
    $ec_num =~ s/-/n/g;
    print OUT_AMB "$ec_num\t$ec_name\n";
  }
}

close (ENZIN);
close (CLSIN);
close (OUT_LIVE);
close (OUT_DEL);
close (OUT_TRANS);
close (OUT_AMB);

