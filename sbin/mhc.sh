#!/bin/bash
ARGML="./foliage"
CPU="select[ model==BL460C_2_26 ]"
BASE="/lustre/scratch1/sanger/lg8/mhc"

A3B7D=$BASE"/A3B7"
A1B8D=$BASE"/A1B8"

FPHD="/fastphase"
LOGD="/log"

A3B7="A3B7"
A1B8="A1B8"

ARGS=2
RS=42
PR=2

submit() {
	#-------------------------------#
	#	$1	job name	#
	#	$2	job dependecy	#
	#	$3	lsf queue	#
	#	$4	log directory	#
	#	$5	command		#
	#-------------------------------#

	DEPEND='';
	[ -z "$2" ] || DEPEND=' -w"'$2'"';

	QUEUE='';
	[ -z "$3" ] || QUEUE=' -q '$3;

	echo 'bsub'$QUEUE' -R"'$CPU'"'$DEPEND' -J"'$1'" -o"'$4$LOGD'/LSF/'$1'.log" -e"'$4$LOGD'/'$1'.log" "'$5'"'|bash
}

submit-array() {
	#-------------------------------#
	#	$1	job name	#
	#	$2	job dependecy	#
	#	$3	lsf queue	#
	#	$4	number of jobs	#
	#	$5	log directory	#
	#	$6	command		#
	#-------------------------------#

	DEPEND='';

	[ -z "$2" ] || DEPEND=' -w"'$2'"';

	QUEUE='';
	[ -z "$3" ] || QUEUE=' -q '$3;

	echo 'bsub'$QUEUE' -R"'$CPU'"'$DEPEND' -J"'$1'[1-'$4']" -o"'$5$LOGD'/LSF/"%I"_'$1'.log" -e"'$5$LOGD'/"%I"_'$1'.log" "'$6'"'|bash
}

#-------------------------------------------------------------------------------#
#	Create log directories (LSF can't seem to create them on its own)	#
# ------------------------------------------------------------------------------#
mkdir $BASE$LOGD;
mkdir $BASE$LOGD'/LSF/';
mkdir $A3B7D$LOGD;
mkdir $A3B7D$LOGD'/LSF/';
mkdir $A1B8D$LOGD;
mkdir $A1B8D$LOGD'/LSF/';

#---------------------------------------#
#	Create haplotype files		#
# --------------------------------------#
A3B7FPHI="MM_7_fastP.inp"
A3B7FPHO="M7_hapguess_switch_CLEAN.out"
A1B8FPHI="MM_8_fastP.inp"
A1B8FPHO="M8_hapguess_switch.out"

submit '37FH' '' 'small' $A3B7D  './foliage from-phase --I '$A3B7D$FPHD' --O '$A3B7D' --fpi '$A3B7FPHI' --fpo '$A3B7FPHO' --o '$A3B7'.mgr'
submit '18FH' '' 'small' $A1B8D  './foliage from-phase --I '$A1B8D$FPHD' --O '$A1B8D' --fpi '$A1B8FPHI' --fpo '$A1B8FPHO' --o '$A1B8'.mgr'

#---------------------------------------#
#	Sub sample the A1B8 set		#
# --------------------------------------#
submit '18HP' 'ended(18FH)' 'small' $A1B8D './foliage subsample-haplotype --I '$A1B8D' --O '$A1B8D' --i '$A1B8'.mgr --sss '$RS' --ss '$PR' --o '$A1B8'.mgr'

#-----------------------------------------------#
#	Infer Ancestral Recombination Graphs	#
# ----------------------------------------------#
submit '37HM' 'ended(37FH)' 'normal' $A3B7D './foliage margarita --I '$A3B7D' --O '$A3B7D' --i '$A3B7'.mgr --ss '$ARGS' | gzip -c > '$A3B7D'/'$A3B7'.arg.gz'
submit '18HM' 'ended(18FH)' 'long' $A1B8D './foliage margarita --I '$A1B8D' --O '$A1B8D' --i '$A1B8'.mgr --ss '$ARGS' | gzip -c > '$A1B8D'/'$A1B8'.arg.gz'
submit-array '18HMP' 'ended(18HP)' 'normal' $PR $A1B8D './foliage margarita --I '$A1B8D' --O '$A1B8D' --i \${LSB_JOBINDEX}_'$A1B8'.mgr --ss '$ARGS' | gzip -c > '$A1B8D'/\${LSB_JOBINDEX}_'$A1B8'.arg.gz'

#-----------------------------------------------#
#	Calculate statistics on the ARGs	#
# ----------------------------------------------#
submit '37S' 'ended(37HM)' 'normal' $A3B7D 'zcat '$A3B7D'/'$A3B7'.arg.gz|./foliage statistic --I '$A3B7D' --O '$A3B7D' --n '$A3B7' --h '$A3B7'.mgr | gzip -c > '$A3B7D'/'$A3B7'.stats.gz'
submit '18S' 'ended(18HM)' 'basement' $A1B8D 'zcat '$A1B8D'/'$A1B8'.arg.gz|./foliage statistic --I '$A1B8D' --O '$A1B8D' --n '$A1B8' --h '$A1B8'.mgr | gzip -c > '$A1B8D'/'$A1B8'.stats.gz'
submit-array '18SP' 'ended(18HMP)' 'normal' $PR $A1B8D 'zcat '$A1B8D'/\${LSB_JOBINDEX}_'$A1B8'.arg.gz|./foliage statistic --I '$A1B8D' --O '$A1B8D' --n \${LSB_JOBINDEX}_'$A1B8' --h \${LSB_JOBINDEX}_'$A1B8'.mgr --o \${LSB_JOBINDEX}_'$A1B8'.stats'

#---------------------------------------------------------------#
#	Merge all A1B8 statistic files into one gz file		#
# --------------------------------------------------------------#
submit '18SPM' 'ended(18SP)' 'normal' $A1B8D './foliage scollect --I '$A1B8D' --O '$A1B8D' --pt [0-9]+_'$A1B8'\.stats | gzip -c > '$A1B8D'/'$A1B8'.subsample.stats.gz'

#-----------------------------------------------#
#	Create summary statistic for A1B8	#
# ----------------------------------------------#
submit '18PSS' 'ended(18SPM)' 'normal' $A1B8D 'zcat '$A1B8D'/'$A1B8'.subsample.stats.gz|./foliage ssummarize --I '$A1B8D' --O '$A1B8D' --n '$A1B8'P --pt [0-9]+_'$A1B8' --o '$A1B8'.subsummary.stats'

#---------------------------------------#
#	Create collection file		#
# --------------------------------------#
submit 'SC1' 'ended(18PSS) && ended(18S)' 'small' $BASE 'zcat '$A1B8D'/'$A1B8'.stats.gz|./foliage sfilter --I '$A1B8D' --O '$A1B8D' --b '$A1B8'.subsummary.stats --o '$A1B8'.both.stats'
submit 'SC2' 'ended(SC1) && ended(37S)' 'small' $BASE 'zcat '$A3B7D'/'$A3B7'.stats.gz|./foliage sfilter --I '$A1B8D' --O '$BASE' --b '$A1B8'.both.stats | gzip -c > '$BASE'/both.stats.gz'

#-----------------------------------------------#
#	tar the files to save space		#
# ----------------------------------------------#
submit '37ZC' 'ended(SC2)' 'small' $BASE './cleanup.sh '$A3B7D' '$A3B7
submit '18ZC' 'ended(SC2)' 'normal' $BASE './cleanup.sh '$A1B8D' '$A1B8

echo 'done'

