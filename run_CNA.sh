#!/bin/sh

#Initial Processing for identifying the existence of Copy number variations (CNVs)
#Depends on: Bowtie, Samtools, GNU Parallel, R package(DNAcopy).
#Author: Min Hu, Darlan Conterno Minussi
#Last Updated: Jul 5th, 2019


main(){
	create_folder
#	run_bowtie_varbin
#	read_align_stat
#	filter
#	segmentation
#	heatmap
	dashboard
#	ratio_plots
	clean
     }

help="
***version 1.4****
Name
        run_CNV.sh -- Call CNV in a population of N samples

Description

        Initial Processing for identifying the existence of Copy number variations (CNVs)

Usage
        sh `basename $0` [--fastq] <fastq_path_file> [--output] <outdir_path> [--sample] <tissue sample name> [--res] <resolution of bin size> [--cpu] < Cpu number> [--undo_prune] <the proportional increase number for the CBS segment function> [alpha] <significance levels for the CBS test [filter_CellWithEmptyBin] < filter the cells with empty reads in bins > [--makeFig] <plot the CNV profile per cell>  1>process.info 2>process.log

arguments:
        -f|--fastq: [ The path file to the fastq ]
	-o|--output: [ output(default) ]
	-s|--sample: [ sample(default) ]
	-r|--res:	[ 100 200(default) 300 400 500 1000 2000 5000 10000 ]
	-c|--cpu: [ 30(default) ]
	-u|--undo_prune: [ 0.05(default) ]
	-a|--alpha: [ 0.0001(default) ]
	-p|--facs: [ facs path ]
	-e|--filter_CellWithEmptyBin: [ 0.1(default) ]
	-m|--makeFig: [ TRUE(default) NO]
Example:
	sh `basename $0` --fastq fastq_input --output output --sample Tn5
or
	sh `basename $0` --fastq fastq_input --output output --sample Tn5 --res 200 --cpu 30 --undo_prune 0.05 --alpha 0.0001 --makeFig "TRUE"
"

if [[ $# -lt 3 ]] ; then echo "$help" ; exit 1 ; fi

TEMP=`getopt -o f:o:s:r:c:u:a:e:m:p: --long fastq:,output:,sample:,res:,cpu:,undo_prune:,alpha:,filter_CellWithEmptyBin:,facs:,makeFig: \
        -n 'example.bash' -- "$@"`
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP"
while true ; do
        case "$1" in
        -f|--fastq)
                fastq_input=$2 ; shift 2 ;;
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
        -p|--facs)
                facs=$2 ; shift 2 ;;
        -e|--filter_CellWithEmptyBin)
                filter_CellWithEmptyBin=$2 ; shift 2 ;;
	-m|--makeFig)
		makeFig=$2 ; shift 2 ;;
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
[ -z "$filter_CellWithEmptyBin" ] && filter_CellWithEmptyBin="0.1"

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
	cat $fastq_input |grep -v "_R2_" | while read line
	do
	r1_input=${line}
	r2_input=`ls $r1_input|sed -e 's/_R1_/_R2_/g'`
	l1_input=`ls $r1_input|sed -e 's/L00[0-9]/L001/g'`
	l2_input=`ls $r1_input|sed -e 's/L00[0-9]/L002/g'`
	l3_input=`ls $r1_input|sed -e 's/L00[0-9]/L003/g'`
	l4_input=`ls $r1_input|sed -e 's/L00[0-9]/L004/g'`
	nu=0
	for fq in $l1_input $l2_input $l3_input $l4_input
	do
	if [[ -f $fq ]];then nu=$((nu+1));fi
	done

	if [[ $nu == 1 ]]
	then
		if [ ! -f $r2_input ] || [ $r1_input = $r2_input ]
		then
		r2_input="";
		fi
		echo "perl $bin/run_bowtie2_varbin.pl -fq1 $r1_input -fq2 $r2_input -samdir $sam_folder -bamdir $bam_folder -sortdir $sort_folder -vb_dir $vbdir_folder -stat_dir $stat_folder -res $res" >> $s1_tmp
	else
	s1_dir=`dirname $r1_input`
	echo "perl $bin/run_bowtie2_varbin_nextseq.pl -fqdir $s1_dir -samdir $sam_folder -bamdir $bam_folder -sortdir $sort_folder -vb_dir $vbdir_folder -stat_dir $stat_folder -res $res " >> $s1_tmp

	fi

	done
	sort $s1_tmp |uniq > $s2_tmp
	gnu_parallel=$(grep "parallel" $lib/CNA.config | cut -d "=" -f 2)
	cpu="$(expr $cpu / 6)"
	$gnu_parallel -j $cpu < $s2_tmp
	rm -f $output/bowtie-*
	rm -f $s1_tmp
	rm -f $s2_tmp
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
        gcinputFile=$lib/$(grep "gc_2M" $lib/CNA.config | cut -d "=" -f 2b|xargs basename)

 elif [[ $res == "5000"  ]];then
        filterFile=$lib/$(grep "excluded_bin_5M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
        gcinputFile=$lib/$(grep "gc_5M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
 else
        filterFile=$lib/$(grep "excluded_bin_10M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
        gcinputFile=$lib/$(grep "gc_10M" $lib/CNA.config | cut -d "=" -f 2|xargs basename)
 fi
        R CMD BATCH --vanilla --args -filter="$filterFile" -gcinputFile="$gcinputFile"  -outdir="$output" -cpu="$cpu" -filter_CellWithEmptyBin="$filter_CellWithEmptyBin" "$bin/filterScript_GC.R" "$logs_folder/filter_log_file.log"
}

segmentation(){
        source1=$bin/$(grep "CNProcessingCentromeres" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
        source2=$bin/$(grep "Mergelevels.R" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
        chrominfo=$lib/$(grep "chrominfo" $lib/CNA.config |cut -d "=" -f 2|xargs basename)
        R CMD BATCH --vanilla --args -sample="$sample" -undoprune="$undo_prune" -alpha="$alpha" -cpu="$cpu" -chrominfo="$chrominfo" -makeFig="$makeFig" -source1="$source1" -source2="$source2" -outdir="$output" "$bin/CN_Seg_GC_filter_correctionRG.R" "$logs_folder/seg_log_file.log"
}

heatmap(){

        R CMD BATCH --args -sample="$sample" -outdir="$output" "$bin/heatmap.R" "$logs_folder/heatmap.log"
}

dashboard(){
	if [[ -e $facs ]]
	then
	Rscript $bin/render_dashboard_facs.R "$root_dir" "$sample" "$output" "$facs"
	else
	Rscript $bin/render_dashboard.R "$root_dir" "$sample" "$output"
	fi

}

ratio_plots(){
  
        Rscript $bin/ratio_plots.R "$bin" "$output"
  
}

clean(){
	find $output/ -name "*.sam" -type f -exec rm {} \;
	rm -rf $output/../.RData $output/../Rplots.pdf
	rm -rf $bam_folder $sam_folder
	find $sort_folder/ -type d |grep -v "./$" |xargs rm -rf
}

 main
