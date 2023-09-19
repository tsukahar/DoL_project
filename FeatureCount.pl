All scripts are available at /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2022-01-10_FeatureCount

The script for Bowtie2-mapped files are pasted below. 



#!/usr/bin/perl

use Time::Local;

use Term::ANSIColor; 





#########################################

$FEATURE_COUNTS="/home/saurabha/UTILITIES/SUBREAD/subread-1.5.0-p1-source/bin/featureCounts";





$FOLDER="/nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2021-12-26_Bowtie2/step3_sort_index";



@FILES=qw(



sorted_DMSO_Bic_1h_rep1.bam:::DMSO_1h_Bic_rep1

sorted_DMSO_Bic_1h_rep2.bam:::DMSO_1h_Bic_rep2

sorted_DMSO_Bic_1h_rep3.bam:::DMSO_1h_Bic_rep3

sorted_DMSO_Bic_4h_rep1.bam:::DMSO_4h_Bic_rep1

sorted_DMSO_Bic_4h_rep2.bam:::DMSO_4h_Bic_rep2

sorted_DMSO_Bic_4h_rep3.bam:::DMSO_4h_Bic_rep3

sorted_DMSO_water_4h_rep1.bam:::DMSO_4h_water_rep1

sorted_DMSO_water_4h_rep2.bam:::DMSO_4h_water_rep2

sorted_DMSO_water_4h_rep3.bam:::DMSO_4h_water_rep3

sorted_MM401_Bic_1h_rep1.bam:::MM401_1h_Bic_rep1

sorted_MM401_Bic_1h_rep2.bam:::MM401_1h_Bic_rep2

sorted_MM401_Bic_1h_rep3.bam:::MM401_1h_Bic_rep3

sorted_MM401_Bic_4h_rep1.bam:::MM401_4h_Bic_rep1

sorted_MM401_Bic_4h_rep2.bam:::MM401_4h_Bic_rep2

sorted_MM401_Bic_4h_rep3.bam:::MM401_4h_Bic_rep3

sorted_MM401_water_1h_rep1.bam:::MM401_1h_water_rep1

sorted_MM401_water_1h_rep2.bam:::MM401_1h_water_rep2

sorted_MM401_water_1h_rep3.bam:::MM401_1h_water_rep3

sorted_MM401_water_4h_rep1.bam:::MM401_4h_water_rep1

sorted_MM401_water_4h_rep2.bam:::MM401_4h_water_rep2

sorted_MM401_water_4h_rep3.bam:::MM401_4h_water_rep3





);

for $x (0 .. $#FILES)

{

($BAMFILE,$MARK)=split(":::",$FILES[$x]);



system ("ln -s $FOLDER/$BAMFILE ./$MARK.bam");

print color ("cyan"),`ls -lh $FOLDER/$BAMFILE`,color ("yellow"),`ls -lh $MARK.bam`,"\n\n";

sleep 2;

$FILELIST=$FILELIST." $MARK.bam";

}

print color ("magenta"),"\n\n\n$FILELIST\n\n\n",color("reset");





$OPTIONS="-T 2 -a /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2022-01-10_FeatureCount/SAF_INPUT_WHOLEGENES_MM10.SAF -F SAF -g Geneid -s 2 -p -B -C";$COMMAND="$FEATURE_COUNTS $OPTIONS $FILELIST -o MERGED_BRUSEQ_WHOLEGENE_MM10_BT2";

system ("$COMMAND");

system("cat MERGED_BRUSEQ_WHOLEGENE_MM10_BT2 |grep -v Program:featureCounts > MERGED_BRUSEQ_WHOLEGENE_MM10_FEATURECOUNTS_OUTPUT_BT2");