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
		"bamdir:s"=>\$bam_dir,
		"sortdir:s"=>\$sortdir,
		"vb200_dir:s"=>\$vb200_dir,
		"stat_dir:s"=>\$stat_dir,

	  );
die `pod2text $0` if ($Help || !$bam_dir || !$sortdir);


my $bin=dirname($0);
my $support=dirname($bin);
my $config="$support\/lib\/CNA.config";
#print "$config\n";
#@files=glob("$fqdir/*_R1*.fastq $fqdir/*_R1*.fastq.gz");
@files=glob("$bam_dir/*.bam");
my $i=0;
my $pre="";
my $bowtie=find_path("$config","bowtie");
my $samtools=find_path("$config","samtools");
my $bowtie_hg19=find_path("$config","bowtie_hg19");
my $varbin_python=find_path("$config","varbin_python");
my $sortdir_sub;
foreach $file (@files)
{
	my $fname=basename($file);
	$fname =~ /^(.*?)\./;
#	$fname =~ /^(.*?)\.(.*?)\.bam/;
	$pre=$1;
	$cmd="$samtools view $bam_dir/$fname >$sortdir/$pre.sort.sam";
	print "$cmd\n";
	system("$cmd");
#Creates varbins files
	$cmd="$varbin_python $sortdir/$pre.sort.sam $vb200_dir\/$pre\.vb200 $stat_dir\/$pre\.200kb\.stat\.txt";
	print "$cmd\n";
	system("$cmd");
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
