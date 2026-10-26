#/!/bin/bash

### GLOBAL PARAMETERS to be exported to all scripts
pwd="/Users/User/Project/proj/260119_CUTandRUN_S2"
prefix="/scratchc/ghlab/gab/proj/260119_CUTandRUN_S2"

digits=3

###############################
# Function that prepares the scripts for any step in the pipeline
# Arguments:
# 	2) Path of script
# 	1) Name of script
###############################

unset done_prepare; declare -A done_prepare

function prepare {
	# Copy each script only once
	if [ "${done_prepare[$2]}" != "1" ]; then
		echo Preparing $2...
		mkdir -p $pwd/$1
		mkdir -p $pwd/$1/log
		cp templates/$2.template $pwd/$1/$2.sh
		sed -i "s|__PREFIX__|$prefix|g" "$pwd/$1/$2.sh"
		sed -i "s|__PWD__|$pwd|g" "$pwd/$1/$2.sh"
		sed -i "s|__DIGITS__|$digits|g" "$pwd/$1/$2.sh"
		# Copy additional scripts and data
		if [ -d "templates/$2" ]; then
			cp -r templates/$2/* $pwd/$1
		fi
		done_prepare[$2]="1"
	fi
}

mkdir -p done										# Finished steps
unset jobids; declare -A jobids				# Technical replicates
unset jobids_te; declare -A jobids_te		# Technical replicates
unset jobids_viruses; declare -A jobids_viruses			# Technical replicates
unset replicates; declare -A replicates	# Biological replicates
unset replicates_te; declare -A replicates_te	# Biological replicates
unset replicates_viruses; declare -A replicates_viruses	# Biological replicates
unset pooled; declare -A pooled				# Pooled replicates
unset macs2; declare -A macs2				# All MACS2 runs

mkdir -p ../analysis
mkdir -p ../analysis/misc
mkdir -p ../QC

################# TECHNICAL REPLICATES ##################
files=`ls ../data/*.fq.gz`

for file in $files
do
   id=`echo $file | cut -d'/' -f3 | cut -d'.' -f1-$((digits+1))`

	# FASTQC
	if [ ! -f "done/01_fastqc_$id.ok" ]; then
		prepare fastqc fastqc
		sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=2000 --kill-on-invalid-dep=yes -t 1:00:00 -J fastqc_$id -o $pwd/fastqc/log/$id.log -e $pwd/fastqc/log/$id.log ../fastqc/fastqc.sh $id
		touch done/01_fastqc_$id.ok
	fi
done

files=`ls ../data/*.fq.gz`

for file in $files
do
   id=`echo $file | cut -d'/' -f3 | cut -d'.' -f1-$((digits))`
	jobid="0"
	jobid2="0"

	# TRIM
	if [ ! -f "done/02_trim_$id.ok" ]; then
		prepare trim trim
		jobid=`sbatch --nodes=1 --cpus-per-task=2 --ntasks=1 -p epyc --mem=4000 --kill-on-invalid-dep=yes -t 3:00:00 -J trim-$id -o $pwd/trim/log/$id -e $pwd/trim/log/$id ../trim/trim.sh $id | cut -d' ' -f4`
		touch done/02_trim_$id.ok
		rm -f done/03_alignment_viruses_$id.ok
		rm -f done/04_alignment_te_$id.ok
      rm -f done/05_alignment_$id.ok
	fi

	# ALIGNMENT VIRUSES - parallel thread
 	if [ ! -f "done/03_alignment_viruses_$id.ok" ]; then
 		prepare alignment_viruses alignment_viruses
 		dep=""; if [ $jobid != 0 ]; then dep="-d afterok:$jobid"; fi
 jobid2=`sbatch $qos --nodes=1 --cpus-per-task=20 --ntasks=1 -p epyc --mem=66000 --kill-on-invalid-dep=yes -J alignv-$id -t 5:00:00 -o $pwd/alignment_viruses/log/$id -e $pwd/alignment_viruses/log/$id $dep ../alignment_viruses/alignment_viruses.sh $id | cut -d' ' -f4`
 		touch done/03_alignment_viruses_$id.ok
 		id2=`echo $id | cut -d'.' -f1-$((digits-1))`
       rm -f done/09_merge_viruses_$id2.ok
 		jobids_viruses[$id]=$jobid2
 	fi

 	# ALIGNMENT TEs - also parallel thread
 	if [ ! -f "done/04_alignment_te_$id.ok" ]; then
 		prepare alignment_te alignment_te
 		dep=""; if [ $jobid != 0 ]; then dep="-d afterok:$jobid"; fi
 		jobid2=`sbatch $qos --nodes=1 --cpus-per-task=20 --ntasks=1 -p epyc --mem=66000 --kill-on-invalid-dep=yes -J alignte-$id -t 5:00:00 -o $pwd/alignment_te/log/$id -e $pwd/alignment_te/log/$id $dep ../alignment_te/alignment_te.sh $id | cut -d' ' -f4`
 		touch done/04_alignment_te_$id.ok
 		id2=`echo $id | cut -d'.' -f1-$((digits-1))`
       rm -f done/10_merge_te_$id2.ok
 		jobids_te[$id]=$jobid2
 	fi

	# ALIGNMENT STEP1
	if [ ! -f "done/05_alignment_$id.ok" ]; then
		prepare alignment alignment
		dep=""; if [ $jobid != 0 ]; then dep="-d afterok:$jobid"; fi
		jobid=`sbatch --nodes=1 --cpus-per-task=20 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 5:00:00 -J align1-$id -o $pwd/alignment/log/$id -e $pwd/alignment/log/$id $dep ../alignment/alignment.sh $id | cut -d' ' -f4`
		touch done/05_alignment_$id.ok
		id2=`echo $id | cut -d'.' -f1-$((digits-1))`
		rm -f done/06_merge_$id2.ok
		jobids[$id]=$jobid
	fi
done

for key in ${!jobids[@]}; do
    echo ${key} ${jobids[${key}]}
done

######## Load private Python environment #########
conda activate deeptools 
######deeptools
##################################################

################ INDIVIDUAL BIOLOGICAL REPLICATES #################
ids=`ls ../data/*.fq.gz | cut -d'/' -f3 | cut -d'.' -f1-$((digits-1)) | sort | uniq`

for id in $ids
do
	jobid=0
	jobid2=0
	jobid3=0

	# MERGE TECHNICAL REPLICATES
	if [ ! -f "done/06_merge_$id.ok" ]; then
		path="../analysis/misc/merge"
		prepare analysis/misc/merge merge
		techs=`ls ../data/$id* | cut -d'/' -f3 | cut -d'.' -f1-$digits`
		unset dep
		for tech in $techs
		do
			if [ -n "${jobids[$tech]}" ]; then
				dep="$dep:${jobids[$tech]}"
			fi
		done
		if [ ! -z ${dep+x} ]; then
			dep="-d afterok$dep"
		fi
		jobid=`sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 3:00:00 -J merge-$id -o $path/log/$id -e $path/log/$id $dep $path/merge.sh $id | cut -d' ' -f4`
		touch done/06_merge_$id.ok
		replicates["$id"]=$jobid	# Save jobids of technical replicates for final jobs
		rm -f done/11_macs2_$id.ok
		rm -f done/07_tracks_$id.ok
		rm -f done/08_phantompeakqualtools_$id.ok
		id2=`echo $id | cut -d'.' -f1-$((digits-2))`
		rm -f done/14_pool_$id2.ok
      rm -f done/11_deeptools_te_$id.ok
		rm -f done/18_alignment_stats.ok
		rm -f done/19_fingerprint.ok
		rm -f done/13_deeptools_ind_$id.ok
	fi

   # MERGE VIRAL ALIGNMENTS
   if [ ! -f "done/09_merge_viruses_$id.ok" ]; then
      path="../analysis/misc/merge_viruses"
      prepare analysis/misc/merge_viruses merge_viruses
      techs=`ls ../data/$id* | cut -d'/' -f3 | cut -d'.' -f1-$digits`
      unset dep
      for tech in $techs
      do
         if [ -n "${jobids_viruses[$tech]}" ]; then
            dep="$dep:${jobids_viruses[$tech]}"
         fi
      done
      if [ ! -z ${dep+x} ]; then
         dep="-d afterok$dep"
      fi
      jobid3=`sbatch $qos --nodes=1 --cpus-per-task=2 --ntasks=1 -p epyc --mem 80000 --kill-on-invalid-dep=yes -J vmerge-$id -t 1:00:00 -o $path/log/$id -e $path/log/$id $dep $path/merge_viruses.sh $id | cut -d' ' -f4`
      touch done/09_merge_viruses_$id.ok
      replicates_viruses["$id"]=$jobid3   # Save jobids of technical replicates for final jobs
      rm -f done/19_countreads_viruses.ok
   fi

   # MERGE TE ALIGNMENTS
   if [ ! -f "done/10_merge_te_$id.ok" ]; then
      path="../analysis/misc/merge_te"
      prepare analysis/misc/merge_te merge_te
      techs=`ls ../data/$id* | cut -d'/' -f3 | cut -d'.' -f1-$digits`
      unset dep
      for tech in $techs
      do
         if [ -n "${jobids_te[$tech]}" ]; then
            dep="$dep:${jobids_te[$tech]}"
         fi
      done
      if [ ! -z ${dep+x} ]; then
         dep="-d afterok$dep"
      fi
      jobid2=`sbatch $qos --nodes=1 --cpus-per-task=2 --ntasks=1 -p epyc --mem 80000 --kill-on-invalid-dep=yes -J temerge-$id -t 1:00:00 -o $path/log/$id -e $path/log/$id $dep $path/merge_te.sh $id | cut -d' ' -f4`
      touch done/10_merge_te_$id.ok
      replicates_te["$id"]=$jobid2  # Save jobids of technical replicates for final jobs
      rm -f done/11_deeptools_te_$id.ok
		rm -f done/20_countreads_te.ok
   fi

	# DEEPTOOLS TE ALIGNMENTS INDIVIDUAL REPLICATES
	if [ ! -f "done/11_deeptools_te_$id.ok" ]; then
		prepare deepTools_te deepTools_te
		path="$pwd/deepTools_te"
		unset dep
		if [ $jobid != 0 ]; then dep="$dep:$jobid"; fi
		if [ $jobid2 != 0 ]; then dep="$dep:$jobid2"; fi
       if [ ! -z ${dep+x} ]; then dep="-d afterok$dep"; fi
		sbatch --nodes=1 --cpus-per-task=4 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 5:00:00 -J dtte-$id -o $path/log/$id -e $path/log/$id $dep $path/deepTools_te.sh $id
		touch done/11_deeptools_te_$id.ok
	fi

done


################ INDIVIDUAL BIOLOGICAL REPLICATES #################
# DO PER NON-INPUT REPLICATE; ALSO CONSIDER INPUT BEING MERGED
#ids=`ls ../data/*.fq.gz | cut -d'/' -f3 | cut -d'.' -f1-$((digits-1)) | sort | uniq | grep -v input`
ids=`ls ../data/*.fq.gz | cut -d'/' -f3 | cut -d'.' -f1-$((digits-1)) | sort | uniq`

for id in $ids
do
	jobid=0

	# REPLICATES MACS2
	if [[ ! -f "done/11_macs2_$id.ok" ]]; then
		prepare macs2 macs2
		path="$pwd/macs2"

		samplename=`echo $id | cut -d'-' -f2-3 | cut -d'.' -f1-$(($digits-1))`
		unset dep
		if [ -n "${replicates[$id]}" ]; then
			dep="$dep:${replicates[$id]}"
		fi
		input="input-$samplename"
		if [[ "$id" == "$input" ]]; then
			first=`echo $input | cut -d'.' -f1`
			last=`echo $input | cut -d'.' -f2`
			new=1
			if [[ $last == 1 ]]; then new=2; fi
			input="$first.$new"
		fi
echo "Sample: $id, Input: $input"
		if [ -n "${replicates[$input]}" ]; then
			dep="$dep:${replicates[$input]}"
		fi
		if [ ! -z ${dep+x} ]; then
			dep="-d afterok$dep"
		fi
		jobid=`sbatch --nodes=1 --cpus-per-task=4 --ntasks=1 -p epyc --mem=80000 -t 4:00:00 --kill-on-invalid-dep=yes -J m2-$id -o $path/log/$id -e $path/log/$id $dep $path/macs2.sh $id | cut -d' ' -f4`
		macs2["$id"]=$jobid	# Save jobids of biological replicates for final jobs
		touch done/11_macs2_$id.ok
#		rm -f done/12_homer_$id.ok
#		rm -f done/12_fold_change.ok
#		rm -f done/22_frip.ok
	fi

#	# HOMER MOTIFS
#	if [ ! -f "done/12_homer_$id.ok" ]; then
#		prepare homer homer
#		path="$pwd/homer"
#		dep=""; if [ $jobid != 0 ]; then dep="-d afterok:$jobid"; fi
#		sbatch --nodes=1 --cpus-per-task=8 --ntasks=1 -p epyc --mem=16000 --kill-on-invalid-dep=yes -J homer-$id -o $path/log/$id -e $path/log/$id $dep $path/homer.sh $id
#		touch done/12_homer_$id.ok
#	fi

	# DEEPTOOLS INDIVIDUAL REPLICATES
	if [ ! -f "done/13_deeptools_ind_$id.ok" ]; then
		prepare deepTools_ind deepTools_ind
		path="$pwd/deepTools_ind"
		unset dep
		if [ -n "${replicates[$id]}" ]; then
			dep="$dep:${replicates[$id]}"
		fi
		if [ ! -z ${dep+x} ]; then
			dep="-d afterok$dep"
		fi
		sbatch --nodes=1 --cpus-per-task=4 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 6:00:00 -J dti-$id -o $path/log/$id -e $path/log/$id $dep $path/deepTools_ind.sh $id
		touch done/13_deeptools_ind_$id.ok
	fi

done


################ POOLED SAMPLES ################

ids=`ls ../data/*.fq.gz | cut -d'/' -f3 | cut -d'.' -f1-$((digits-2)) | sort | uniq`
for id in $ids
do
	jobid=0

	# POOL BIOLOGICAL REPLICATES
	if [ ! -f "done/14_pool_$id.ok" ]; then
		path="../analysis/misc/pool"
		prepare analysis/misc/pool pool
		bios=`ls ../data/$id* | cut -d'/' -f3 | cut -d'.' -f1-$((digits-1))`
		unset dep
		for bio in $bios
		do
			if [ -n "${replicates[$bio]}" ]; then
				dep="$dep:${replicates[$bio]}"
			fi
		done
		if [ ! -z ${dep+x} ]; then
			dep="-d afterok$dep"
		fi
		jobid=`sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 2:00:00 -J pool-$id -o $path/log/$id -e $path/log/$id $dep $path/pool.sh $id | cut -d' ' -f4`
		touch done/14_pool_$id.ok
		rm -f done/15_tracks_pooled_$id.ok
		rm -f done/17_macs2_pooled_$id.ok
		rm -f done/16_deeptools_$id.ok
		pooled["$id"]=$jobid
	fi

done

################ POOLED SAMPLES ################
# DO PER NON-INPUT POOL, ALSO CONSIDER INPUT BEING MERGED
ids=`ls ../data/*.fq.gz | cut -d'/' -f3 | cut -d'.' -f1-$((digits-2)) | sort | uniq`

for id in $ids
do
	jobid=0

#	# POOLED MACS2
#	if [[ ! -f "done/17_macs2_pooled_$id.ok" ]]; then
#		prepare macs2_pooled macs2_pooled
#		path="$pwd/macs2_pooled"
#
#		samplename=`echo $id | cut -d'-' -f2-3 | cut -d'.' -f1-$(($digits-2))`
#		unset dep
#		if [ -n "${pooled[$id]}" ]; then
#			dep="$dep:${pooled[$id]}"
#		fi
#		input="input-$samplename"
#echo "Sample: $id, Input: $input"
#		if [ -n "${pooled[$input]}" ]; then
#			dep="$dep:${pooled[$input]}"
#		fi
#		if [ ! -z ${dep+x} ]; then
#			dep="-d afterok$dep"
#		fi
#
#		jobid=`sbatch --nodes=1 --cpus-per-task=4 --ntasks=1 -p epyc --mem=8000 --kill-on-invalid-dep=yes -J m2p-$id -o $path/log/$id -e $path/log/$id $dep $path/macs2_pooled.sh $id | cut -d' ' -f4`
#		macs2["$id"]=$jobid	# Save jobids of technical replicates for final jobs
#		touch done/17_macs2_pooled_$id.ok
#		rm -f done/16_deeptools_$id.ok
#		rm -f done/21_fold_changes.ok
#		rm -f done/22_frip.ok
#	fi

	# DEEPTOOLS
	if [ ! -f "done/16_deeptools_$id.ok" ]; then
		prepare deepTools deepTools
		path="$pwd/deepTools"
		unset dep
		if [ -n "${pooled[$id]}" ]; then
			dep="$dep:${pooled[$id]}"
		fi
		if [ ! -z ${dep+x} ]; then
			dep="-d afterok$dep"
		fi
		sbatch --nodes=1 --cpus-per-task=4 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 6:00:00 -J dtp-$id -o $path/log/$id -e $path/log/$id $dep $path/deepTools.sh $id
		touch done/16_deeptools_$id.ok
	fi
done

######## Unload private Python environment ########
conda deactivate
###################################################

################ ALL TECHNICAL REPLICATES #################

# Define a list of all replicates
unset all_techs
for key in ${!jobids[@]}; do
    all_techs="$all_techs:${jobids[${key}]}"
done
if [ ! -z ${all_techs+x} ]; then
	all_techs="-d afterok$all_techs"
fi

# SUMMARY OF ALIGNMENT STATS
if [ ! -f "done/18_alignment_stats.ok" ]; then
	path="../analysis/misc/alignment_stats"
	prepare analysis/misc/alignment_stats alignment_stats
	sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=80000 --kill-on-invalid-dep=yes -t 1:00:00 -J stats -o $path/log/log.txt -e $path/log/log.txt $all_techs $path/alignment_stats.sh
	touch done/18_alignment_stats.ok
fi

################ ALL BIOLOGICAL REPLICATES ##################

# Define a list of all VIRUS replicates
unset all_replicates
for key in ${!replicates_viruses[@]}; do
    all_replicates="$all_replicates:${replicates_viruses[${key}]}"
done
if [ ! -z ${all_replicates+x} ]; then
   all_replicates="-d afterok$all_replicates"
fi

# COUNTREADS VIRUSES
if [ ! -f "done/19_countreads_viruses.ok" ]; then
   prepare countreads_viruses countreads_viruses
   path="../countreads_viruses"
   jobid=`sbatch $qos --nodes=1 --cpus-per-task=20 --mem 80000 --ntasks=1 -p epyc -t 2:00:00 --kill-on-invalid-dep=yes -J count -o $path/log/log.txt -e $path/log/log.txt $all_replicates $path/countreads_viruses.sh $id | cut -d' ' -f4`
   touch done/19_countreads_viruses.ok
fi

# Define a list of all TE replicates
unset all_replicates
for key in ${!replicates_te[@]}; do
    all_replicates="$all_replicates:${replicates_te[${key}]}"
done
if [ ! -z ${all_replicates+x} ]; then
   all_replicates="-d afterok$all_replicates"
fi

# COUNTREADS TE
if [ ! -f "done/20_countreads_te.ok" ]; then
   prepare countreads_te countreads_te
   path="../countreads_te"
   jobid=`sbatch $qos --nodes=1 --cpus-per-task=20 --mem 80000 --ntasks=1 -p epyc -t 2:00:00 --kill-on-invalid-dep=yes -J count -o $path/log/log.txt -e $path/log/log.txt $all_replicates $path/countreads_te.sh $id | cut -d' ' -f4`
   touch done/20_countreads_te.ok
fi

# Define a list of all replicates
unset all_replicates
for key in ${!replicates[@]}; do
    all_replicates="$all_replicates:${replicates[${key}]}"
done
if [ ! -z ${all_replicates+x} ]; then
	all_replicates="-d afterok$all_replicates"
fi

conda activate deeptools
# SUMMARY OF ALIGNMENT STATS
if [ ! -f "done/19_fingerprint.ok" ]; then
	path="../analysis/misc/fingerprint"
	prepare analysis/misc/fingerprint fingerprint
	sbatch --nodes=1 --cpus-per-task=20 --ntasks=1 -p epyc --mem=64G --kill-on-invalid-dep=yes -t 5:00:00 -J fp -o $path/log/log.txt -e $path/log/log.txt $all_replicates $path/fingerprint.sh
	touch done/19_fingerprint.ok
fi
conda deactivate

################ ALL MACS2 RUNS ##################
# Define a list of all replicates
unset all_macs2
for key in ${!macs2[@]}; do
    all_macs2="$all_macs2:${macs2[${key}]}"
done
if [ ! -z ${all_macs2+x} ]; then
	all_macs2="-d afterok$all_macs2"
fi

# MACS2 FOLD CHANGES STATS
#if [ ! -f "done/21_fold_changes.ok" ]; then
#	path="../QC/fold_changes"
#	prepare QC/fold_changes fold_changes
#	sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=2000 --kill-on-invalid-dep=yes -J m2_FC -o $path/log/log.txt -e $path/log/log.txt $all_macs2 $path/fold_changes.sh
#	touch done/21_fold_changes.ok
#fi

# MACS2 FRIP
#if [ ! -f "done/22_frip.ok" ]; then
#	path="../QC/FRIP"
#	prepare QC/FRIP frip
#	sbatch --nodes=1 --cpus-per-task=1 --ntasks=1 -p epyc --mem=2000 --kill-on-invalid-dep=yes -J frip -o $path/log/log.txt -e $path/log/log.txt $all_macs2 $path/frip.sh
#	touch done/22_frip.ok
#fi


