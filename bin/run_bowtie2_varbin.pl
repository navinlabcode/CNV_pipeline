#! /usr/bin/perl -w
use File::Basename;
use Cwd;
use Getopt::Long;


=head1 Name

d1 Name


=head1 Description

generate sample alignment based on bowtie alignment


=head1 Usage

perl run_bowtie2.pl [-fqdir] <fq_path> [-outdir]<output>
--verbose   output running progress information to screen
--help      output help information to screen

=head1 Example

perl run_bowtie2.pl -fqdir fq_path -outdir output -bedfile bed_name

=cut

###get options from command line into variables and set default values
my ($Verbose,$Help,$fqdir,$outdir,$bedfile);
GetOptions(
		"verbose"=>\$Verbose,
		"help"=>\$Help,
		"fq1:s"=>\$fq1,
		"fq2:s"=>\$fq2,
		"samdir:s"=>\$samdir,
		"bamdir:s"=>\$bamdir,
		"sortdir:s"=>\$sortdir,
		"vb_dir:s"=>\$vb_dir,
		"stat_dir:s"=>\$stat_dir,
		"res:s"=>\$res
	  );
die `pod2text $0` if ($Help || !$fq1 || !$bamdir || !$samdir);


my $bin=dirname($0);
my $support=dirname($bin);
my $config="$support\/lib\/CNA.config";

my $i=0;

my $bowtie=find_path("$config","bowtie");
my $samtools=find_path("$config","samtools");
my $bowtie_hg19=find_path("$config","bowtie_hg19");
my $varbin_python=find_path("$config","varbin_python");
$varbin_python=$bin."/".`basename $varbin_python`;
chomp $varbin_python;

my $chrominfo=find_path("$config","chrominfo");
my $bins;
if ($res =~/^200$/){
	$bins=find_path("$config","bins_200k");
}elsif($res=~/^100$/){
	$bins=find_path("$config","bins_100k");
}elsif($res=~/^300$/){
	$bins=find_path("$config","bins_300k");
}elsif($res=~/^400$/){
	$bins=find_path("$config","bins_400k");
}elsif($res=~/^500$/){
	$bins=find_path("$config","bins_500k");
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


my $sortName2;my $sortName;



my $fname=basename($fq1);
#$fname =~ /^(.*?)\_(.*?)\_(.*?)\_(.*?)\_(.*?)\.fastq/;
my $pre;
if($fq2){
	$fname =~ /^(.*?)\_R\d_001/;        #MR1-MDAMB231-C100_S484_L002_R1_001.fastq.gz
		$pre=$1;
}else{
	$fname=~/^(.*?)\.fastq/;
	$pre=$1;
}
$pre=~s/\_L00\d//g;
$samName=$samdir."/".$pre.".sam";

#get sam file
my $cmd_sam;
if($fq2){
	$cmd_sam="$bowtie -x $bowtie_hg19/hg19 -1 $fq1 -2 $fq2 -S $samName -p 6 ";
}else{
	$cmd_sam="gzip -dc $fq1 |$bowtie -x $bowtie_hg19/hg19 - -S $samName -p 6 ";
}

#from sam to bam
$bamName=$bamdir."/".$pre.".bam";
$cmd_bam="$samtools view -bS -q 1 $samName > $bamName";

$sortName=$sortdir."/".$pre.".sort";
$cmd_sort="$samtools sort  $bamName $sortName";
$sortName2="$sortName.bam";
system("$cmd_sam");
system("$cmd_bam");
system("$cmd_sort");


#Converts sorted Bam file to a sorted Sam file using Samtools 
$cmd="$samtools view $sortName2 > $sortName.sam";
system("$cmd");

#Creates varbins files
$cmd="$varbin_python $sortName\.sam $vb_dir\/$pre\.vb $stat_dir\/$pre\.stat\.txt $chrominfo $bins";
print "$cmd\n";
system("$cmd");


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
