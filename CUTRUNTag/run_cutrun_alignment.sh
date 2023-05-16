############################################################################
##
##  Usage:
##      Cutrun_pipeline.sh configure-file
##
##  Juliana Lee (juliana_lee@hms.harvard.edu)
##  CAUTION: NOT intended for F1 use
##
############################################################################

#!/bin/bash

if [ $# -ne 1 ]; then
        echo "$0 <Configure file>"
    exit
fi

CONFIG=$1


###############################################################################
### Set Enviroments
###############################################################################


if [ -f $CONFIG ]
then
    sed -e 's/ *$//' $CONFIG > tmp.txt
    mv tmp.txt $CONFIG
    dos2unix $CONFIG

    FLIST=`grep -w '^FLIST' $CONFIG | cut -d '=' -f2`
    FASTQ_DIR=`grep -w '^FASTQ_DIR' $CONFIG | cut -d '=' -f2`
    LOG_DIR=`grep -w '^LOG_DIR' $CONFIG | cut -d '=' -f2`
    OUTPUT_DIR=`grep -w '^OUTPUT_DIR' $CONFIG | cut -d '=' -f2`
    SCRIPT_DIR=`grep -w '^SCRIPT_DIR' $CONFIG | cut -d '=' -f2`
    FASTQ_SUFFIX=`grep -w '^FASTQ_SUFFIX' $CONFIG | cut -d '=' -f2`

    GENOME_INDEX_B6=`grep -w '^GENOME_INDEX_B6' $CONFIG | cut -d '=' -f2`
    GENOME_INDEX_CAST=`grep -w '^GENOME_INDEX_CAST' $CONFIG | cut -d '=' -f2`

    Bowtie2_PARAMS=`grep -w '^Bowtie2_PARAMS' $CONFIG | cut -d '=' -f2`

    MIN_QUALITY_SCORE=`grep -w '^MIN_QUALITY_SCORE' $CONFIG | cut -d '=' -f2`

    MOD_B6=`grep -w '^MOD_B6' $CONFIG | cut -d '=' -f2`
    MOD_CAST=`grep -w '^MOD_CAST' $CONFIG | cut -d '=' -f2`

    BLACKLIST=`grep -w '^BLACKLIST' $CONFIG | cut -d '=' -f2`
    GRCm38_GENE_ANNOTATION=`grep -w '^GRCm38_GENE_ANNOTATION' $CONFIG | cut -d '=' -f2`
    GRCm38_effective_genome_size=`grep -w '^GRCm38_effective_genome_size' $CONFIG | cut -d '=' -f2`

    GCC=`grep -w '^GCC' $CONFIG | cut -d '=' -f2`
    GSL=`grep -w '^GSL' $CONFIG | cut -d '=' -f2`
    PYTHON=`grep -w '^PYTHON' $CONFIG | cut -d '=' -f2`
    JAVA=`grep -w '^JAVA' $CONFIG | cut -d '=' -f2`
    Conda=`grep -w '^Conda' $CONFIG | cut -d '=' -f2`
    Samtools=`grep -w '^Samtools' $CONFIG | cut -d '=' -f2`
    Picard=`grep -w '^Picard' $CONFIG | cut -d '=' -f2`
    Deeptools=`grep -w '^Deeptools' $CONFIG | cut -d '=' -f2`
    Betools=`grep -w '^Betools' $CONFIG | cut -d '=' -f2`
    RSeQC=`grep -w '^RSeQC' $CONFIG | cut -d '=' -f2`
    Multiqc=`grep -w '^Multiqc' $CONFIG | cut -d '=' -f2`

    CORE=`grep -w '^CORE' $CONFIG | cut -d '=' -f2`
    CORE_Bowtie=`grep -w '^CORE_Bowtie' $CONFIG | cut -d '=' -f2`
    CORE_Samtools=`grep -w '^CORE_Samtools' $CONFIG | cut -d '=' -f2`
    CORE_lapels=`grep -w '^CORE_lapels' $CONFIG | cut -d '=' -f2`

else
    echo "Please provide a configuration file"
    exit
fi



mkdir -p $LOG_DIR
mkdir -p $OUTPUT_DIR


###############################################################################
### Pipeline steps
###############################################################################
ARRAY_JOBID=()

while IFS=$'\t' read -r -a Array
do
       sample_name=${Array[0]}
       sample_ID=${Array[1]}
       batch=${Array[2]}
       owner=${Array[3]// /_}
       species=${Array[4]}

       echo $sample_name


##############################################################################
# Merge technical replicates/lanes
# ############################################################################
cat $FASTQ_DIR/${sample_ID}_*_R1_*.fastq.gz > $FASTQ_DIR/${sample_ID}.R1.fastq.gz
cat $FASTQ_DIR/${sample_ID}_*_R2_*.fastq.gz > $FASTQ_DIR/${sample_ID}.R2.fastq.gz

FASTQS="$FASTQ_DIR/${sample_ID}.R1.fastq.gz $FASTQ_DIR/${sample_ID}.R2.fastq.gz"

mkdir -p $OUTPUT_DIR/$sample_ID

#############################################################################
# Adapter and quality trimming
#############################################################################
rm -f $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
cp $SCRIPT_DIR/sbatch.sh $OUTPUT_DIR/$sample_ID/adapter_trimming.sh

echo -e "#SBATCH --mem-per-cpu=4G" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "#SBATCH -c $CORE" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "#SBATCH -o $LOG_DIR/${sample_ID}.adapter_trimming.log" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "module load gcc/6.2.0 python/3.6.0 $Conda" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e 'eval "$(conda shell.bash hook)"' >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "conda activate bowtie_conda" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "module load trimgalore/0.6.6" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "trim_galore --cores $CORE -o $OUTPUT_DIR/$sample_ID --fastqc --paired $FASTQS" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh
echo -e "conda deactivate" >> $OUTPUT_DIR/$sample_ID/adapter_trimming.sh


#############################################################################
# Read Alignment
#############################################################################
rm -f $OUTPUT_DIR/$sample_ID/alignment.sh
cp $SCRIPT_DIR/sbatch.sh $OUTPUT_DIR/$sample_ID/alignment.sh

echo -e "#SBATCH --mem-per-cpu=20G" >> $OUTPUT_DIR/$sample_ID/alignment.sh
echo -e "#SBATCH -c $CORE_Bowtie" >> $OUTPUT_DIR/$sample_ID/alignment.sh
echo -e "#SBATCH -o $LOG_DIR/${sample_ID}.alignment.log" >> $OUTPUT_DIR/$sample_ID/alignment.sh
echo -e "module load gcc/6.2.0 python/3.6.0 bowtie2 $Conda" >> $OUTPUT_DIR/$sample_ID/alignment.sh

echo -e "bowtie2 $Bowtie2_PARAMS --threads $CORE_Bowtie -x $GENOME_INDEX_B6 -1 $OUTPUT_DIR/$sample_ID/${sample_ID}.R1_val_1.fq.gz -2 $OUTPUT_DIR/$sample_ID/${sample_ID}.R2_val_2.fq.gz -S $OUTPUT_DIR/$sample_ID/${sample_ID}.sam" >> $OUTPUT_DIR/$sample_ID/alignment.sh
    
#############################################################################
# Post Alignment Filtering
#############################################################################
rm -f $OUTPUT_DIR/$sample_ID/samtools.sh
cp $SCRIPT_DIR/sbatch.sh $OUTPUT_DIR/$sample_ID/samtools.sh

echo -e "#SBATCH --mem-per-cpu=16G" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "#SBATCH -c $CORE_Samtools" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "#SBATCH -o $LOG_DIR/${sample_ID}.samtools.log" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "module load $GCC" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "module load $Samtools" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "module load $JAVA" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "module load $Betools" >> $OUTPUT_DIR/$sample_ID/samtools.sh

echo -e "samtools view -b $OUTPUT_DIR/$sample_ID/${sample_ID}.sam > $OUTPUT_DIR/$sample_ID/${sample_ID}.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "samtools view -bh -f 3 -F 12 -q $MIN_QUALITY_SCORE $OUTPUT_DIR/$sample_ID/${sample_ID}.bam > $OUTPUT_DIR/$sample_ID/${sample_ID}.mapped.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "samtools view -h $OUTPUT_DIR/$sample_ID/${sample_ID}.mapped.bam | sed '/chrM/d;/random/d;/chrUn/d' | samtools sort -@ $CORE_Samtools -O bam -o $OUTPUT_DIR/$sample_ID/${sample_ID}.filtered.sorted.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "java -Xms32000m -jar /n/groups/cbdm-db/ly82/ATACSeq_pipeline/tools/picard.jar MarkDuplicates I=$OUTPUT_DIR/$sample_ID/${sample_ID}.filtered.sorted.bam O=$OUTPUT_DIR/$sample_ID/${sample_ID}.filtered.sorted.nodup.bam M=$OUTPUT_DIR/$sample_ID/${sample_ID}.dup.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "bedtools intersect -abam $OUTPUT_DIR/$sample_ID/${sample_ID}.filtered.sorted.nodup.bam -b $BLACKLIST -wa -v > $OUTPUT_DIR/$sample_ID/${sample_ID}.final.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "cp $OUTPUT_DIR/$sample_ID/${sample_ID}.filtered.sorted.nodup.bam $OUTPUT_DIR/$sample_ID/${sample_ID}.final.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "samtools sort $OUTPUT_DIR/$sample_ID/${sample_ID}.final.bam > $OUTPUT_DIR/$sample_ID/${sample_ID}.final.sorted.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "samtools index $OUTPUT_DIR/$sample_ID/${sample_ID}.final.sorted.bam" >> $OUTPUT_DIR/$sample_ID/samtools.sh
echo -e "rm -f $OUTPUT_DIR/$sample_ID/${sample_ID}.sam" >> $OUTPUT_DIR/$sample_ID/samtools.sh

##################################################
# bedgraph and bigwig
##################################################
rm -f $OUTPUT_DIR/$sample_ID/convert_format.sh
cp $SCRIPT_DIR/sbatch.sh $OUTPUT_DIR/$sample_ID/convert_format.sh

echo -e "#SBATCH --mem-per-cpu=4G" >> $OUTPUT_DIR/$sample_ID/convert_format.sh
echo -e "#SBATCH -c $CORE" >> $OUTPUT_DIR/$sample_ID/convert_format.sh
echo -e "#SBATCH -o $LOG_DIR/${sample_ID}.convert_format.log" >> $OUTPUT_DIR/$sample_ID/convert_format.sh
echo -e "module load $GCC" >> $OUTPUT_DIR/$sample_ID/convert_format.sh
echo -e "module load $PYTHON" >> $OUTPUT_DIR/$sample_ID/convert_format.sh
echo -e "module load $Deeptools" >> $OUTPUT_DIR/$sample_ID/convert_format.sh

echo -e "bamCoverage --bam $OUTPUT_DIR/$sample_ID/${sample_ID}.final.sorted.bam --outFileName $OUTPUT_DIR/$sample_ID/${sample_ID}.bw --outFileFormat bigwig --binSize 1 --normalizeUsing CPM --effectiveGenomeSize $GRCm38_effective_genome_size --numberOfProcessors $CORE" >> $OUTPUT_DIR/$sample_ID/convert_format.sh

echo -e "bamCoverage --bam $OUTPUT_DIR/$sample_ID/${sample_ID}.final.sorted.bam --outFileName $OUTPUT_DIR/$sample_ID/${sample_ID}.bdg --outFileFormat bedgraph --binSize 1 --normalizeUsing CPM --effectiveGenomeSize $GRCm38_effective_genome_size --numberOfProcessors $CORE" >> $OUTPUT_DIR/$sample_ID/convert_format.sh

ARRAY_JOBID+=("$CONVERT_FORMAT_JOBID")

done < $FLIST

cd $OUTPUT_DIR

for i in $OUTPUT_DIR/*; do
    cd $i
    sbatch adapter_trimming.sh
    ADAPTER_TRIMMING_ID=$(sbatch --parsable adapter_trimming.sh)
    sbatch --parsable --dependency=afterok:$ADAPTER_TRIMMING_ID alignment.sh
    ALIGNMENT_JOBID=$(sbatch --parsable --dependency=afterok:$ADAPTER_TRIMMING_ID alignment.sh)
    sbatch --parsable --dependency=afterok:$ALIGNMENT_JOBID samtools.sh
    SAMTOOLS_JOBID=$(sbatch --parsable --dependency=afterok:$ALIGNMENT_JOBID samtools.sh)
    sbatch --parsable --dependency=afterok:$SAMTOOLS_JOBID convert_format.sh
    cd .. ;
done



