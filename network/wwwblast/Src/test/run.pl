#! /usr/local/bin/perl

if (not -e "Log") {
  `mkdir Log`;
}

my $app = $ARGV[0];
my $diff = 'diff -w';
my $time = &GetTimeCmd();

#my $oldbin = "/net/blast012/export/home/web/public/htdocs/BLAST/bl2seq/$app";
my $oldbin = "../wblast2.old.REAL";
chomp(my $basedir = `pwd`);
my $newbin = "$basedir/$app";

my $out = "out";

if (not -e "$out") {
  `mkdir $out`;
}

my %Tests;

if ($app eq "wblast2.REAL") {
   $Tests{'blastp'} = "\"ONE=129295&TWO=XP_222492.2&FILTER=1&PROGRAM=blastp\"";

   $Tests{'blastn'} = "\"ONE=555&TWO=101&FILTER=1&PROGRAM=blastn\"";

   $Tests{'megablast'} = "\"ONE=555&TWO=AC091728&FILTER=1&PROGRAM=blastn&MEGABLAST=yes&WORD=20\"";

   $Tests{'tblastn'} = "\"ONE=9930103&TWO=9930102&FILTER=1&PROGRAM=tblastn\"";
 
   $Tests{'blastx'} = "\"ONE=3090&TWO=3091&FILTER=1&PROGRAM=blastx\"";

   $Tests{'tblastx'} = "\"ONE=555&TWO=101&PROGRAM=tblastx&FILTER=1&WORD=3\"";

   $Tests{'blastn-minus'} = "\"ONE=NT_004487.15&TWO=AA441981.1&FROM=7685545&TO=7686027&FFROM=10&TTO=480&STRAND=2&FILTER=1&PROGRAM=blastn\"";

   $Tests{'blastn-plus'} = "\"ONE=NT_004487.15&TWO=AA441981.1&FROM=7685545&TO=7686027&FFROM=10&TTO=480&STRAND=1&FILTER=1&PROGRAM=blastn\"";

   $Tests{'fully-masked'} = "\"ONE=U09816&TWO=BX641126.1&FROM=1280&TO=1324&FFROM=2052&TTO=2082&STRAND=2&PROGRAM=blastn\"";
} else {
  if ($app eq "blast_cs.REAL") {


  }
}
foreach $test (keys %Tests) {
    print "\nTest ", $test, "";
    print "\n----------------\n";
    print "Parameters: $Tests{$test}\n\n";

    foreach $binary_type (qw(New Old)) {

        print "\t\"$binary_type\".  Time: ";

        if ($binary_type eq "New") {
            $binary = $newbin;
        } else {
            $binary = $oldbin;
        }
        $rv = system("$time sh -c '$binary $Tests{$test}' > $out/$test.$binary_type.out 2> $out/$test.$binary_type.err");
        $time_str = `tail -3 $out/$test.$binary_type.err | tr -s "\n" " "`;
        chomp($time_str);
        print $time_str, "\n";
    }
    `$diff $out/$test.Old.out $out/$test.New.out > $out/$test.diff`;
}

sub GetTimeCmd() {
    return "/usr/bin/time -p" if (`uname` =~ /linux/i);
    return "/usr/bin/time";
}

