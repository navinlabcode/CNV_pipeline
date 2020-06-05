#! /usr/bin/perl -w
use File::Basename;
use Cwd;
use Getopt::Long;


=head1 Name

d1 Name

run_bowtie2_pe_SeqCap_EZV2.pl  -- bowtie alignment,mark the duplication reads and caculate coverage

=head1 Description

generate sample alignment based on bowtie alignment


=head1 Usage

perl run_bowtie2.pl [-fqdir] <fq_path> [-outdir]<output>
--verbose   output running progress information to screen
--help      output help information to screen

=head1 Exmple

perl run_bowtie2.pl -fqdir fq_path -outdir output -bedfile bed_name

=cut

###get options from command line into variables and set default values
my ($Verbose,$Help,$fqdir,$outdir,$bedfile);
GetOptions(
		"verbose"=>\$Verbose,
		"help"=>\$Help,
		"samdir:s"=>\$sam_dir,
		"vb_dir:s"=>\$vb_dir,
		"stat_dir:s"=>\$stat_dir,
		"res:s"=>\$res,

	  );
die `pod2text $0` if ($Help || !$sam_dir);


my $bin=dirname($0);
my $support=dirname($bin);
my $config="$support\/lib\/CNA.config";
my $i=0;
my $pre="";
my $bowtie=find_path("$config","bowtie");
my $samtools=find_path("$config","samtools");
my $bowtie_hg19=find_path("$config","bowtie_hg19");
my $chrominfo=find_path("$config","chrominfo");
my $varbin_python=find_path("$config","varbin_python");
$varbin_python=$bin."/".`basename $varbin_python`;
chomp $varbin_python;
my $bins;

if ($res =~/^200$/){
 $bins=find_path("$config","bins_200k");
}elsif($res=~/^1000$/){
 $bins=find_path("$config","bins_1M");
}elsif($res=~/^2000$/){
 $bins=find_path("$config","bins_2M");
}elsif($res=~/^5000$/){
 $bins=find_path("$config","bins_5M");
}elsif($res=~/^10000$/){
 $bins=find_path("$config","bins_10M");
}

$bins=$support."/lib/".`basename $bins`;
chomp $bins;
$chrominfo=$support."/lib/".`basename $chrominfo`;
chomp $chrominfo;

my $sortdir_sub;
open FIL,$sam_dir or die $!;
my @files=<FIL>;
foreach $file (@files)
{
	chomp $file;
	my $fname=basename($file);
	$fname =~ /^(.*?)\./;
#	$fname =~ /^(.*?)\.(.*?)\.bam/;
	$pre=$1;
#	$cmd="$samtools view $sam_dir/$fname >$sortdir/$pre.sort.sam";
#	print "$cmd\n";
#	system("$cmd");
#Creates varbins files
	$cmd="$varbin_python $file $vb_dir\/$pre\.vb $stat_dir\/$pre\.stat\.txt $chrominfo $bins";
	print "$cmd\n";
	#system("$cmd");
}


sub find_path{
	my($config,$software_name)=@_;
	open CON, $config or die $!;
	while (<CON>){
		chomp;
		next unless(/\=/);
		my($name,$path)=split(/\=/,$_);
		if($name=~/$software_name/){
			return $path;		
		}	
	}
	close CON;
}
