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

my $input="oncoarray_bcac_public_release_oct17.chr$chr.txt";
open F2, $input;
while(<F2>){
my @line=split/\s+/;
my($chr,$pos,$ref,$alt)=split/_/,$line[0];
$ref=uc $ref;
$alt=uc $alt;
my $id="$chr-$pos-$ref-$alt";
if(defined $data{$id}){
$line[1]=$data{$id};
$line[4]=uc $line[4];
$line[5]=uc $line[5];
}
my $print= join("\t",@line);
print "$print\n";
}

