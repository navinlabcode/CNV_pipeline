#!/bin/sh

#Initial Processing for identifying the existence of Copy number variations (CNVs)
#Depends on: Bowtie, Samtools, GNU Parallel, R package(DNAcopy).
#Author: Min Hu, Darlan Conterno Minussi
#Last Updated: Oct 15th, 2021


main(){
	create_folder
	run_bowtie_varbin
	read_align_stat
	filter
	segmentation
	ratio_plots
	dashboard
	clean
}

help="
***version 1.4****
Name
run_CNV.sh -- Call CNV in a population of N samples

Description

Initial Processing for identifying the existence of Copy number variations (CNVs)

	Usage
	sh `basename $0` [--fastq] <fastq_path_file> [--sam] <sam_file> [--output] <outdir_path> [--sample] <tissue sample name> [--res] <resolution of bin size> [--cpu] < Cpu number> [--undo_prune] <the proportional increase number for the CBS segment function> [alpha] <significance levels for the CBS test [filter_CellWithEmptyBin] < filter the cells with empty reads in bins > [--makeFig] <plot the CNV profile per cell>  1>process.info 2>process.log

	arguments:
	-f|--fastq: [ The path file to the fastq ]
	-l|--sam: [sam input]
	-b|--bam: [bam input]
	-o|--output: [ output(default) ]
	-s|--sample: [ sample(default) ]
	-r|--res:	[ 100 200(default) 300 400 500 1000 2000 5000 10000 ]
	-c|--cpu: [ 30(default) ]
	-u|--undo_prune: [ 0.05(default) ]
	-a|--alpha: [ 0.0001(default) ]
	-p|--facs: [ facs path ]
	-e|--filter_CellWithEmptyBin: [ 0.1(default) ]
	-t|--filter_ReadCount: [ 100000(default) ]
	-m|--makeFig: [ TRUE(default) NO]
	-g|--genome_len: [ ref genome length file path]
	-q|--miseq : [ False(default) TRUE]
	Example:
	sh `basename $0` --fastq fastq_input --output output --sample Tn5
	or
	sh `basename $0` --fastq fastq_input --filter_ReadCount 100 --filter_CellWithEmptyBin 0.9 --miseq True [if your data yield is very low]
	or
	sh `basename $0` --sam sam_input --genome_len /volumes/seq/code/PIPELINES/CNA_pipeline_v1.4/lib/genome.fa.fai
	or
	sh `basename $0` --bam bam_input --output output --sample Tn5
	or
	sh `basename $0` --fastq fastq_input --output output --sample Tn5 --res 200 --cpu 30 --undo_prune 0.05 --alpha 0.0001 --makeFig "TRUE"
	"

	if [[ $# -lt 3 ]] ; then echo "$help" ; exit 1 ; fi

	TEMP=`getopt -o f:l:b:o:s:r:c:u:a:p:e:t:m:g:q: --long fastq:,sam:,bam:,output:,sample:,res:,cpu:,undo_prune:,alpha:,facs:,filter_CellWithEmptyBin:,filter_ReadCount:,makeFig:,genome_len:,miseq: \
       -n 'example.bash' -- "$@"`
       if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
       eval set -- "$TEMP"
       while true ; do
       case "$1" in
       -f|--fastq)
       fastq_input=$2 ; shift 2 ;;
       -l|--sam)
       sam_input=$2 ; shift 2 ;;
       -b|--bam)
       bam_input=$2 ; shift 2 ;;
       -o|--output)
       output=$2 ; shift 2 ;;
       -s|--sample)
       sample=$2 ; shift 2 ;;
       -r|--res)
       res=$2 ; shift 2 ;;
       -c|--cpu)
       cpu=$2 ; shift 2 ;;
       -u|--undo_prune)
       undo_prune=$2 ; shift 2 ;;
       -a|--alpha)
       alpha=$2 ; shift 2 ;;
       -g|--genome_len)
	genome_len=$2 ; shift 2 ;;
       -p|--facs)
       facs=$2 ; shift 2 ;;
       -e|--filter_CellWithEmptyBin)
       filter_CellWithEmptyBin=$2 ; shift 2 ;;
       -t|--filter_ReadCount)
       filter_ReadCount=$2 ; shift 2 ;;
       -m|--makeFig)
       makeFig=$2 ; shift 2 ;;
       -q|--miseq)
       miseq=$2 ; shift 2 ;;
       --) shift ; break ;;
       *) echo "Internal error!" ; exit 1 ;;
       esac
       done

       [ -z "$output" ] && output="output"
       [ -z "$sample" ] && sample="sample"
       [ -z "$res" ] && res="200"
       [ -z "$cpu" ] && cpu="30"
       [ -z "$undo_prune" ] && undo_prune="0.05"
       [ -z "$alpha" ] && alpha="0.0001"
       [ -z "$makeFig" ] && makeFig="TRUE"
       [ -z "$miseq" ] && miseq="FALSE"
       [ -z "$filter_CellWithEmptyBin" ] && filter_CellWithEmptyBin="0.1"
       [ -z "$filter_ReadCount" ] && filter_ReadCount="100000"
       [ -f $fastq_input ] || { echo " Warning !!! $fastq_input does not exist Please check your fastq files !!! " ;echo "$help";exit 1; }


create_folder(){
	root_dir=`dirname $0`
	dashboard=`dirname $0`/dashboard
	bin=`dirname $0`/bin
	lib=`dirname $0`/lib
	sam_folder="$output/sam"
	bam_folder="$output/bam"
	sort_folder="$output/sort"
	vbdir_folder="$output/vbdir"
	stat_folder="$output/stat"
	metrics_folder="$output/metrics"
	logs_folder="$output/logs"
	final_folder="$output/final_result"
	ratio_plots_folder="$output/final_result/ratio_plots"
	mkdir -p $sam_folder/
	mkdir -p $bam_folder/
	mkdir -p $sort_folder/
	mkdir -p $vbdir_folder/
	mkdir -p $metrics_folder/
	mkdir -p $stat_folder/
	mkdir -p $logs_folder/
	mkdir -p $ratio_plots_folder/
	mkdir -p $final_folder/
    }

run_bowtie_varbin(){
	[ -e "$output/split_bowtie.sh" ] && `rm -f "$output/split_bowtie.sh"`
		s1_tmp="$output/split_bowtie.sh"
		s2_tmp="$output/split_bowtie_uniq.sh"

		if [[ -f $sam_input ]];then
       		[ -z "$genome_len" ] && genome_len=$(grep "genome_len" $lib/CNA.config | cut -d "=" -f 2)
		perl $bin/run_bowtie2_varbin_sam.pl  -sam $sam_input -bamdir $bam_folder -sortdir $sort_folder -vb_dir $vbdir_folder -stat_dir $stat_folder -res $res -genome_len $genome_len >>$s1_tmp
		elif [[ -f $bam_input ]];then 
			perl $bin/run_bowtie2_varbin_bam.pl  -bamdir $bam_input -vb_dir $vbdir_folder -stat_dir $stat_folder -res $res -sortdir $sort_folder >>$s1_tmp
		else
		cat $fastq_input |grep -v "_R2_" | while read line
		do
			r1_input=${line}
		r2_input=`ls $r1_input|sed -e 's/_R1_/_R2_/g'`
			s1_dir=`dirname $r1_input`
			echo "perl $bin/run_bowtie2_varbin.pl -fqdir $s1_dir -samdir $sam_folder -bamdir $bam_folder -sortdir $sort_folder -vb_dir $vbdir_folder -stat_dir $stat_folder -res $res " >> $s1_tmp
			done
		fi
			sort $s1_tmp |uniq > $s2_tmp
			gnu_parallel=$(grep "parallel" $lib/CNA.config | cut -d "=" -f 2)
			cpu="$(expr $cpu / 6)"
			$gnu_parallel -j $cpu < $s2_tmp
#			rm -f $s1_tmp
#			rm -f $s2_tmp
			time=`date`
			echo "$time step1 run_bowtie & varbin  is done"
}


read_align_stat(){
	perl $bin/All_stat_metrics.pl $stat_folder > "$metrics_folder/all_stat_metrics.txt"
	perl $bin/All_stat_metrics_summary.pl $stat_folder > "$metrics_folder/all_stat_metrics_summary.txt"
	perl $bin/All_raw_bin_metrics.pl $vbdir_folder > "$metrics_folder/all_raw_bins_metrics.txt"
}



filter(){
	if [[ $res == "100" ]];then
	filterFile=$lib/$(grep "excluded_bin_100k" $lib/CNA.config | cut -d "=" -f 2 |xargs basename)
	gcinputFile=$lib/$(grep "gc_100k" $lib/CNA.config | cut -d "=" -f 2 |xargs basename)
	elif [[ $res == "200" ]];then
	filterFile=$lib/$(grep "excluded_bin_200k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_200k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "300" ]];then
	filterFile=$lib/$(grep "excluded_bin_300k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_300k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "400" ]];then
	filterFile=$lib/$(grep "excluded_bin_400k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_400k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "500" ]];then
	filterFile=$lib/$(grep "excluded_bin_500k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_500k" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "1000"  ]];then
	filterFile=$lib/$(grep "excluded_bin_1M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_1M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "2000"  ]];then
	filterFile=$lib/$(grep "excluded_bin_2M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_2M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	elif [[ $res == "5000"  ]];then
	filterFile=$lib/$(grep "excluded_bin_5M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_5M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	else
	filterFile=$lib/$(grep "excluded_bin_10M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	gcinputFile=$lib/$(grep "gc_10M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
	fi
	R CMD BATCH --vanilla --args -filter="$filterFile" -gcinputFile="$gcinputFile"  -outdir="$output" -cpu="$cpu" -filter_ReadCount="$filter_ReadCount" -filter_CellWithEmptyBin="$filter_CellWithEmptyBin" "$bin/filterScript_GC.R" "$logs_folder/filter_log_file.log"
}

segmentation(){
	source1=$bin/$(grep "CNProcessingCentromeres" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
		source2=$bin/$(grep "Mergelevels.R" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
		chrominfo=$lib/$(grep "chrominfo" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
		R CMD BATCH --vanilla --args -sample="$sample" -undoprune="$undo_prune" -alpha="$alpha" -cpu="$cpu" -chrominfo="$chrominfo" -makeFig="$makeFig" -source1="$source1" -source2="$source2" -outdir="$output" "$bin/CN_Seg_GC_filter_correctionRG.R" "$logs_folder/seg_log_file.log"
}


dashboard(){
	if [[ -e $facs ]];then
		Rscript $bin/render_dashboard_facs.R "$root_dir" "$sample" "$output" "$facs"
	elif [[ $miseq == "TRUE" ]];then
		Rscript $bin/render_dashboard_miseq.R "$root_dir" "$sample" "$output"
		echo "miseq run"
		echo $miseq 
	else
		Rscript $bin/render_dashboard.R "$root_dir" "$sample" "$output"
	fi

}

ratio_plots(){
	Rscript $bin/ratio_plots.R "$bin" "$output"
}

clean(){
        find $output/ -name "*.sam" -type f |xargs rm -rf
        find $output/ -name "*.sort.bam" -type f |xargs rm -rf
        find $output/ -name "bam" -type d |xargs rm -rf
        find $output/ -name "sam" -type d |xargs rm -rf
        find $output/sort  -type d |awk 'NR>1'| xargs rm -rf
}

main
