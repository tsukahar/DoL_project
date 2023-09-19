#align fastq files with Bowtie2                        

bowtie2 -p 12 -q --local \
-x /nfs/value/siwase2/GENOME_INDICES/MM10/Bowtie2/MM10 \
-1 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-7_CG$
-2 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-7_CG$
-S /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2021-12-26_Bowtie2/MM401_4h_w$

bowtie2 -p 12 -q --local \
-x /nfs/value/siwase2/GENOME_INDICES/MM10/Bowtie2/MM10 \
-1 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-14_G$
-2 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-14_G$
-S /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2021-12-26_Bowtie2/MM401_4h_w$

bowtie2 -p 12 -q --local \
-x /nfs/value/siwase2/GENOME_INDICES/MM10/Bowtie2/MM10 \
-1 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-21_C$
-2 /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2020_input_fastq/2002-MP-21_C$
-S /nfs/dataden/umms-siwase/SHIGEKI/Takao_Bru-seq/2021-12-26_Bowtie2/MM401_4h_w$


