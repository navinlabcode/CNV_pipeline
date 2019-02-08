#!/usr/bin/perl

use strict;

my$input_folder=$ARGV[0];
my@files=glob("$input_folder/*.vb");

my(%hasha,%pos,@pos,@bin_g,@samples_id);
print "chr\tposition\tabs_position\t";
foreach my$file(@files){
	my$sample_name;
	if($file=~/vb(.+)?\/(.+)?\.vb/){
		$sample_name=$2;
		open FIL,"$file" or die $!;
		my @bin_counts;
		while(<FIL>){
			chomp;
			my@array=split(/\s+/,$_);
			my $bin_pos="$array[0]\t$array[1]\t$array[2]";
			if(!exists $pos{$bin_pos}){
				push @pos,$bin_pos;
				$pos{$bin_pos}++;
			}
			
			push @{$hasha{$sample_name}},$array[3];
		}
			push @samples_id,$sample_name;
		close FIL;
		push @bin_g,$hasha{$sample_name};
	}
}

my $all_sample_id=join("\t",@samples_id);
print "$all_sample_id\n";

my$i=0;
foreach my$bin(@pos){
	print "$bin";
	foreach my$sid(@samples_id){
	print "\t@{$hasha{$sid}}[$i]";
	}
	print "\n";
	$i++;
}
