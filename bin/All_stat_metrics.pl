#!/usr/bin/perl

use strict;

my$input_folder=$ARGV[0];
my@files=glob("$input_folder/*.stat.txt");

my(@total_reads,@dups,@readskept,@Medianbincount);

print "Sample Name\tTotalReads\tDupsRemoved\tReadsKept\tMedianBinCount\n";
foreach my$file(@files){
	my$sample_name;
	if($file=~/stat\/(.+)?\.stat/){
		$sample_name=$1;
	}elsif($file=~/([^\.\/]+)\_S\d+\_L001\./){
		$sample_name=$1;
	}
	open FIL,"$file" or die $!;
	while(<FIL>){
		chomp;
		next if(/Total/);
		my@array=split(/\s+/,$_);
#		next if($array[0]<100000);
		push @total_reads,$array[0];
		push @dups,$array[1];
		push @readskept,$array[2];
		push @Medianbincount,$array[3];
		print "$sample_name\t$array[0]\t$array[1]\t$array[2]\t$array[3]\n";
	}
	close FIL;
}

#print "\n\n";
#print "TotalReads,TotalReads Coefficient Of Variation,Mean_TotalReads,Mean_DupsRemoved,Mean_ReadsKept,Mean_BinCount,Mean_Dupreads_Fraction(%)\n";

my$mean_tr=int(&average(@total_reads));
my$mean_du=int(&average(@dups));
my$mean_rk=int(&average(@readskept));
my$mean_bc=int(&average(@Medianbincount));
my$mean_du_fq=($mean_du/$mean_tr);
$mean_du_fq=sprintf "%.4f",$mean_du_fq;
$mean_du_fq*=100;

my $standard_deviation=&stdev(@total_reads);
my $CV=($standard_deviation/$mean_tr);
$CV=sprintf "%.2f",$CV;
my $total_rd=int(&sum(@total_reads));
#print "$total_rd,$CV,$mean_tr,$mean_du,$mean_rk,$mean_bc,$mean_du_fq\%\n";

sub sum{	
        my@data = @_;
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $sum=$total;
        return $sum;
}

sub average{
        my@data = @_;
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $average = $total / @data;
        return $average;
}

sub stdev{
        my@data = @_;
        if($#data == 1){
                return 0;
        }
        my $average = &average(@data);
        my $sqtotal = 0;
        foreach(@data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@data-1)) ** 0.5;
        return $std;
}
