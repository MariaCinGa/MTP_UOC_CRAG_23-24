#!/usr/local/bin/perl
use Getopt::Long;
use File::Find;
use strict;

my (%options,@GENO1,@GENO2,@TABLE1,@TABLE2, @vector, $head1, $head2, @headt,$head_names);
my ($x, $i, $n, $len1, $len2);

#options
GetOptions(
    'i1=s'  => \$options{'i1'},
    'i2=s'  => \$options{'i2'},
    't1=s'  => \$options{'t1'},
    't2=s'  => \$options{'t2'},
    'n1=s'  => \$options{'n1'},
    'n2=s'  => \$options{'n2'},
    'o1=s'  => \$options{'o1'},
    'o2=s'  => \$options{'o2'}
);

#Usage
if($options{'i1'} eq undef || $options{'i2'} eq undef || $options{'t1'} eq undef  || $options{'t2'} eq undef  || $options{'n1'} eq undef || $options{'n2'} eq undef || $options{'o1'} eq undef || $options{'o2'} eq undef) {
   print "\nUsage: \nperl merge_geno_table.pl -i1 (geno_file_pop1)  -i2 (geno_file_pop2) -t1 (table_snps_file_pop1) -t2 (table_snps_file_pop2) -n1 (number_ind1) -n2 (number_ind2) -o1 (output file geno) -o2 (output file table)\n\n";
   exit;
}

#read the geno file 1
$len1=0;
open(INPUT1,"$options{'i1'}") or die;
$_=<INPUT1>;
#$head1=$_; #keep header
while(<INPUT1>) {
    $_ =~ /(\S+)\s(\d+)\s(.+)/;
    if($_ ne "\n") {
        @vector = ($1,$2,$3);
        if($len1>0 && $vector[1]==$GENO1[($len1-1)*3+1]) { #erase mhits
            splice @GENO1,-3;
            $len1--;
        }
        else {
            push (@GENO1,@vector); #keep all in the array GENO1
            $len1++;
        }
    }
}
close(INPUT1);

#read the geno file 2
$len2=0;
open(INPUT2,"$options{'i2'}") or die;
$_=<INPUT2>;
#$_ =~ /(\S+)\s(\S+)\s(.+)/; #keep header (except Chr and Pos)
#$head2=($3);
while(<INPUT2>) {
    $_ =~ /(\S+)\s(\d+)\s(.+)/;
    if($_ ne "\n") {
        @vector = ($1,$2,$3);
        if($len2>0 && $vector[1]==$GENO2[($len2-1)*3+1]) { #erase mhits
            splice @GENO2,-3;
            $len2--;
        }
        else {
            push (@GENO2,@vector); #keep all in the array GENO2
            $len2++;
        }
    }
}
close(INPUT2);

#read the table file 1
open(INPUT1,"$options{'t1'}") or die;
for($i=0;$i<27;$i++) {$_=<INPUT1>; push(@headt,$_)} #keep 27 rows + header
$_=<INPUT1>; #keep header
$head_names=$_;
$head_names =~ s/\R//g;
while(<INPUT1>) {
    $_ =~ /(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/;#\s(.+)/;
    @vector = ($1,$2,$3,$4,$5,$6,$7);
    push (@TABLE1,@vector); #keep all in the array @TABLE1
}
close(INPUT1);

#read the table file 2
open(INPUT2,"$options{'t2'}") or die;
for($i=0;$i<28;$i++) {$_=<INPUT2>;} #skip 27 + header rows
while(<INPUT2>) {
    $_ =~ /(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/;#\s(.+)/;
    @vector = ($1,$2,$3,$4,$5,$6,$7);
    push (@TABLE2,@vector); #keep all in the array @TABLE2
}
close(INPUT2);

#Read SNP POS in pop 1 and 2 and combine both together
#write geno combined file
open(OUTPUT1,">$options{'o1'}") or die;
$head1="Chr Pos";
for($n=1;$n<=$options{'n1'};$n++) {$head1=$head1." P1.IND$n";}
$head2="";
for($n=1;$n<=$options{'n1'};$n++) {$head2=$head2." P2.IND$n";}
print OUTPUT1 $head1.$head2."\n";
$i=0;
for($x=0; ($x<$len1); $x++) {
    while($GENO1[$x*3+1]>$GENO2[$i*3+1] && $i<$len2) {
        print OUTPUT1 $GENO2[$i*3+0]." ".$GENO2[$i*3+1]." ";
        for($n=0;$n<$options{'n1'};$n++) {
            print OUTPUT1 "0 ";
        }
        print OUTPUT1 $GENO2[$i*3+2]."\n";
        $i++;
    }
    if($GENO1[$x*3+1]==$GENO2[$i*3+1]) {
        print OUTPUT1 $GENO1[$x*3+0]." ".$GENO1[$x*3+1]." ".$GENO1[$x*3+2].$GENO2[$i*3+2]."\n";
        $i++;
    }
    else {
        print OUTPUT1 $GENO1[$x*3+0]." ".$GENO1[$x*3+1]." ",$GENO1[$x*3+2];
        for($n=0;$n<$options{'n2'};$n++) {
            print OUTPUT1 "0 ";
        }
        print OUTPUT1 "\n";
    }
}
while($i<$len2) {
    print OUTPUT1 $GENO2[$i*3+0]." ".$GENO2[$i*3+1]." ";
    for($n=0;$n<$options{'n1'};$n++) {
        print OUTPUT1 "0 ";
    }
    print OUTPUT1 $GENO2[$i*3+2]."\n";
    $i++;
}
close(OUTPUT1);

#write table combined file
open(OUTPUT2,">$options{'o2'}") or die;
for($i=0;$i<27;$i++) {print OUTPUT2 $headt[$i];}
print OUTPUT2 $head_names."P1"."\t"."FREQP2\n";
$i=0;
for($x=0; ($x<($#TABLE1+1)/7); $x++) {
    while($TABLE1[$x*7+0]>$TABLE2[$i*7+0] && ($i<($#TABLE2+1)/7)) {
        print OUTPUT2 $TABLE2[$i*7+0]."\t".$TABLE2[$i*7+1]."\t".$TABLE2[$i*7+2]."\t".$TABLE2[$i*7+3]."\t".$TABLE2[$i*7+4]."\t".$TABLE2[$i*7+5]."\t"."0.00"."\t".$TABLE2[$i*7+6]."\n";
        $i++;
    }
    if($i<($#TABLE2+1)/7 && $TABLE1[$x*7+0]==$TABLE2[$i*7+0]) {
        print OUTPUT2 $TABLE1[$x*7+0]."\t".$TABLE1[$x*7+1]."\t".$TABLE2[$i*7+2]."\t".$TABLE1[$x*7+3]."\t".$TABLE1[$x*7+4]."\t".$TABLE1[$x*7+5]."\t".$TABLE1[$x*7+6]."\t".$TABLE2[$i*7+6]."\n";
        $i++;
    }
    else {
        print OUTPUT2 $TABLE1[$x*7+0]."\t".$TABLE1[$x*7+1]."\t".$TABLE1[$x*7+2]."\t".$TABLE1[$x*7+3]."\t".$TABLE1[$x*7+4]."\t".$TABLE1[$x*7+5]."\t".$TABLE1[$x*7+6]."\t"."0.00"."\n";
    }
}
while($i<($#TABLE2+1)/7) {
    print OUTPUT2 $TABLE2[$i*7+0]."\t".$TABLE2[$i*7+1]."\t".$TABLE2[$i*7+2]."\t".$TABLE2[$i*7+3]."\t".$TABLE2[$i*7+4]."\t".$TABLE2[$i*7+5]."\t"."0.00"."\t".$TABLE2[$i*7+6]."\n";
    $i++;
}
close(OUTPUT2);
