#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

sub usage{
	print STDERR "\e[32;1m
Description
	Calling variant from a single mpileup file and output vcf.

Usage
	perl $0 -i file.mpileup -o file.vcf -s sample_name -t threshuold_for_snp

Parameters
	-i <s> input file in mpileup format
	-o <s> output file in vcf-like format
	-s <s> sample name used in vcf header line
	-t <s> threshold for snp, about e-3 for illumina reads.[default 0.005]

Author
	Bingke Jiao bkjiao\@genetics.ac.cn
Date
	2023.7
	\e[0m\n";
}

my ($input,$output,$sample,$thre);
GetOptions(
	"i:s" => \$input,
	"o:s" => \$output,
	"s:s" => \$sample,
	"t:s" => \$thre
);

unless($input && $output && $sample){
	&usage;
	exit 1;
}
$thre ||= 0.005;

open OUT,">",$output or die "Cannot open $output\n";
print OUT "#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample\t$sample\_indel_freq\t$sample\_snp_freq\n";

open IN,"<",$input or die "Cannot open $input\n";
while(<IN>){
	chomp;
	my @tmp=split /\t/,$_;
	my ($DP,$bases)=@tmp[3,4];
	my (%hash_site,%hash_indel,%hash_snp);
	my @array=split //,$bases; #Store the variants into an array.
	while($bases=~m/([0-9]+)/g){
		my $pos=pos($bases);
		$pos=$pos-length($1)-1;
		my $indel=substr($bases,$pos,1+length($1)+$1);
		$indel=uc($indel);
		$hash_site{$indel}++;
		$hash_indel{$indel}++;
		for(my $i=0;$i<1+length($1)+$1;$i++){
			$array[$pos+$i]=""; #Assign the null value to the indels.
		}
	
	}
	$bases=join("",@array)."\n";
	#print "The \$bases is $bases\n";
	while($bases=~/([ATCGN])/ig){
		my $snp=uc($1);
		$hash_site{$snp}++;
		$hash_snp{$snp}++;
	}
	my $ref_num=0;
	while($bases=~/[\.,]/g){
		$ref_num++;
	}

	my $sign=1;
	my $snp_num=0;
	foreach my $k (keys %hash_snp){
		$snp_num+=$hash_snp{$k};
	}

	my $indel_num=0;
	foreach my $indel (keys %hash_indel){
		if($hash_indel{$indel}>0){
			$indel_num+=$hash_indel{$indel};
		}
	}

	if($indel_num > 0){
		$sign=1;
	}else{
		if($DP >0 && $snp_num/$DP >= $thre){
			$sign=1;
		}
	}
	
	
	if($sign==1){
		print OUT join("\t",@tmp[0..2]);
		my @array_value;
		print OUT "\t".join(",",sort(keys %hash_site));
		print OUT "\t.\t.\tDP=$DP\tAD\t";
		foreach my $k (sort keys %hash_site){
			push @array_value,$hash_site{$k};
		}
		print OUT join(",",$ref_num,@array_value);
		if($DP>0){
			printf OUT "\t%.7f\t%.7f\n",$indel_num/$DP,$snp_num/$DP;
		}else{
			printf OUT "\t%.7f\t%.7f\n",0,0;
		}
	}
}
close IN;
close OUT;
