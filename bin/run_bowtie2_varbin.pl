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
		"fqdir:s"=>\$fqdir,
		"samdir:s"=>\$samdir,
		"bamdir:s"=>\$bamdir,
		"sortdir:s"=>\$sortdir,
		"vb_dir:s"=>\$vb_dir,
		"stat_dir:s"=>\$stat_dir,
		"res:s"=>\$res
	  );
die `pod2text $0` if ($Help || !$fqdir || !$bamdir || !$samdir);


my $bin=dirname($0);
my $support=dirname($bin);
my $config="$support\/lib\/CNA.config";
#print "$config\n";
@files=glob("$fqdir/*_R1*.fastq $fqdir/*_R1*.fastq.gz");


my $i=0;
my $pre="";
my $sub_pre="";

my $bowtie=find_path("$config","bowtie");
my $samtools=find_path("$config","samtools");
my $bowtie_hg19=find_path("$config","bowtie_hg19");
my $varbin_python=find_path("$config","varbin_python");
$varbin_python=$bin."/".`basename $varbin_python`;
chomp $varbin_python;
my $chrominfo=find_path("$config","chrominfo"); 
my $sortdir_sub;

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

foreach $file (@files)
{
	my $fname=basename($file);
#$fname =~ /^(.*?)\_(.*?)\_(.*?)\_(.*?)\_(.*?)\.fastq/;
	my $fq1=$file;
	my $fq2=$file;
	$fq2=~s/\_R1\_/\_R2\_/;
	if(-e $fq2){
	$fname =~ /^(.*?)\_R\d_001/;
	$pre=$1;  
	}else{
	$fname =~ /^(.*?)\.fastq/;
	$pre=$1;
	}
	$samName=$samdir."/".$pre.".sam";

#get sam file
	if(-e $fq2){
	$cmd_sam="$bowtie -x $bowtie_hg19/hg19 -1 $fq1 -2 $fq2 -S $samName -p 6 "; 
	}else{
	$cmd_sam="gzip -dc $file |$bowtie -x $bowtie_hg19/hg19 - -S $samName -p 6 ";
	}
#from sam to bam
	$bamName=$bamdir."/".$pre.".bam";
	$cmd_bam="$samtools view -bS -q 1 $samName > $bamName";

#sort bam
	$sub_pre=$pre;
	$sub_pre=~s/\_L00\d//g;
	$sortdir_sub=$sortdir."/".$sub_pre;
	`mkdir -p $sortdir_sub`;
	$sortName=$sortdir_sub."/".$pre.".sorted";
	$cmd_sort="$samtools sort  $bamName $sortName";
#$sortName2="$sortName.bam";
	print STDOUT "$cmd_sam\n";
	system("$cmd_sam");
	print STDOUT "$cmd_bam\n";
	system("$cmd_bam");
	print STDOUT "$cmd_sort\n";
	system("$cmd_sort");
}


@files=glob("$sortdir_sub/*.sorted.bam");

$size=@files;
if($size>1 )
{
	$cmd="$samtools merge $sortdir/$sub_pre.sort.bam ";

	foreach $file (@files)
	{
		$cmd=$cmd." ".$file;
	}

	print STDOUT "$cmd\n";
	print STDOUT "more than one lane in this sample $cmd\n";
	system("$cmd");

}
elsif($size == 1 )
{
	$cmd="mv $files[0] $sortdir/$sub_pre.sort.bam ";
	print STDOUT "$cmd\n";
	system("$cmd");
}
#Converts sorted Bam file to a sorted Sam file using Samtools 
$cmd="$samtools view $sortdir/$sub_pre.sort.bam >$sortdir/$sub_pre.sort.sam";
system("$cmd");
#Creates varbins files
$cmd="$varbin_python $sortdir/$sub_pre.sort.sam $vb_dir\/$sub_pre\.vb $stat_dir\/$sub_pre\.stat\.txt $chrominfo $bins";
print STDOUT "$cmd\n";
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
