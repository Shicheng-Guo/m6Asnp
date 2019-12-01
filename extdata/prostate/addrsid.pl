use strict;
use Cwd;
my $chr=shift @ARGV;
my $dbsnp153="/gpfs/home/guosa/hpc/db/dbSNP153/dbSNP153.hg19.chr$chr.vcf";
my %data;
open F1,$dbsnp153 || die "cannot open $dbsnp153";
while(<F1>){
next if /^#/;
my ($chr,$pos,$rs,$ref,$alt)=split/\s+/;
$ref=uc $ref;
$alt= uc $alt;
my $id="$chr-$pos-$ref-$alt";
$data{$id}=$rs;
}

my $input="meta_v3_onco_euro_overall_ChrAll_1_release.chr$chr.txt";
open F2, $input;
while(<F2>){
print $_ if /MarkerName/;
my @line=split/\s+/;
my($chr,$pos,$ref,$alt)=split/_/,$line[0];
$ref=uc $ref;
$alt=uc $alt;
my $id="$chr-$pos-$ref-$alt";
if(defined $data{$id}){
$line[2]=$data{$id};
$line[5]=uc $line[5];
$line[6]=uc $line[6];
}
my $print= join("\t",@line);
print "$print\n";
}
