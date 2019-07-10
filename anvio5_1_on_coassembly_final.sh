#!/bin/bash
#SBATCH -c 18                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=200gb                         # mem in gb
#SBATCH -t 10-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J anvio5_1.180                  	# Name of job
#SBATCH -e /global/scratch/hpc3565/anvio5_v13_v10.err  
#SBATCH -o /global/scratch/hpc3565/anvio5_v13_v10.out


# processing back on CAC
source activate anvio5

# making environmental variables
R1s=`ls /global/scratch/hpc3565/coassembly/*1P*`
R2s=`ls /global/scratch/hpc3565/coassembly/*2P*`

BASEDIR=/global/scratch/hpc3565
#these contigs have only prokaryotic reads filtered by EukRep
OUTPUT=$BASEDIR/output_anvio5_1
#DNA=$HOME/rna_dna_2017/dna/raw_data/data_complete/output_anvio
#CONTIGS=$DNA/contigs.fa
CONTIGS=$BASEDIR/megahit_coassembly_v13/prokarya_all.fa
READS=$HOME/rna_dna_2017/dna/raw_data/data_complete
MAPPING=$OUTPUT/mapping
KAIJU=$HOME/rna_dna_2017/dna/raw_data/data_complete/kaiju/bin
KAIJUFILES=$HOME/rna_dna_2017/dna/raw_data/data_complete/kaijudb_prokaryotes

######
echo "#######"
echo "Message from Graham: This is run using anvio5, with a sse3 requirement trying to beat the queue" 
echo "#######"
#######

cd $OUTPUT
#simplifying contig file ### using the EukRep output file that should have removed eukaryotes
anvi-script-reformat-fasta $CONTIGS -o $OUTPUT/contigs-fixed.fa -l 1000 --simplify-names --report-file $OUTPUT/name_conversions.txt
mv $OUTPUT/contigs-fixed.fa $OUTPUT/contigs.fa

#building index   
bowtie2-build $OUTPUT/contigs.fa $MAPPING/contigs --threads 18
#indexed BAM file for each samples
#Bowtie2 build 
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_AS3_R_1P.fastq.gz -2 $READS/HI.DNA_AS3_R_2P.fastq.gz -S $MAPPING/HI.DNA_AS3.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_BS3_R_1P.fastq.gz -2 $READS/HI.DNA_BS3_R_2P.fastq.gz -S $MAPPING/HI.DNA_BS3.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_RS3_R_1P.fastq.gz -2 $READS/HI.DNA_RS3_R_2P.fastq.gz -S $MAPPING/HI.DNA_RS3.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_S1_R_1P.fastq.gz -2 $READS/HI.DNA_S1_R_2P.fastq.gz -S $MAPPING/HI.DNA_S1.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_S2_R_1P.fastq.gz -2 $READS/HI.DNA_S2_R_2P.fastq.gz -S $MAPPING/HI.DNA_S2.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_S3_R_1P.fastq.gz -2 $READS/HI.DNA_S3_R_2P.fastq.gz -S $MAPPING/HI.DNA_S3.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_S4_R_1P.fastq.gz -2 $READS/HI.DNA_S4_R_2P.fastq.gz -S $MAPPING/HI.DNA_S4.sam
bowtie2 --threads 18 -x $MAPPING/contigs -1 $READS/HI.DNA_S5_R_1P.fastq.gz -2 $READS/HI.DNA_S5_R_2P.fastq.gz -S $MAPPING/HI.DNA_S5.sam
#samtool conversion sam to bam
samtools view -F 4 -bS $MAPPING/HI.DNA_AS3.sam > $MAPPING/HI.DNA_AS3-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_BS3.sam > $MAPPING/HI.DNA_BS3-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_RS3.sam > $MAPPING/HI.DNA_RS3-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_S1.sam > $MAPPING/HI.DNA_S1-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_S2.sam > $MAPPING/HI.DNA_S2-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_S3.sam > $MAPPING/HI.DNA_S3-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_S4.sam > $MAPPING/HI.DNA_S4-RAW.bam
samtools view -F 4 -bS $MAPPING/HI.DNA_S5.sam > $MAPPING/HI.DNA_S5-RAW.bam
#reformating bam for anvio
anvi-init-bam $MAPPING/HI.DNA_AS3-RAW.bam -o $MAPPING/HI.DNA_AS3.bam
anvi-init-bam $MAPPING/HI.DNA_BS3-RAW.bam -o $MAPPING/HI.DNA_BS3.bam
anvi-init-bam $MAPPING/HI.DNA_RS3-RAW.bam -o $MAPPING/HI.DNA_RS3.bam
anvi-init-bam $MAPPING/HI.DNA_S1-RAW.bam -o $MAPPING/HI.DNA_S1.bam
anvi-init-bam $MAPPING/HI.DNA_S2-RAW.bam -o $MAPPING/HI.DNA_S2.bam
anvi-init-bam $MAPPING/HI.DNA_S3-RAW.bam -o $MAPPING/HI.DNA_S3.bam
anvi-init-bam $MAPPING/HI.DNA_S4-RAW.bam -o $MAPPING/HI.DNA_S4.bam
anvi-init-bam $MAPPING/HI.DNA_S5-RAW.bam -o $MAPPING/HI.DNA_S5.bam
#removing old sam and bam files
rm $MAPPING/HI.DNA_*.sam $MAPPING/HI.DNA_*-RAW.bam

###########
#formulate the database
###########
anvi-gen-contigs-database -f $OUTPUT/contigs.fa -o $OUTPUT/contigs.db -n 'hazen 2017 prokaryote DNA coassembled contigs db'

#hmm profiles
anvi-run-hmms -c $OUTPUT/contigs.db --num-threads 18

#view contigs stats locally by transferring over contigs.db file
anvi-display-contigs-stats $OUTPUT/contigs.db --report-as-text -o $OUTPUT/contig_stats.txt

# annonate genes with NCBI COGs, look into importing other functions 
# if never run before uncomment the line below
#anvi-setup-ncbi-cogs
anvi-run-ncbi-cogs -c $OUTPUT/contigs.db --num-threads 18

######
# echo "have things run fine up to here without an error" 
echo "#######"
echo "Message from Graham: have things run fine up to here without an error?" 
echo "#######"
#######

#create profiles for each sample
#below didn't work so will adjust to do it manually
#for sample in `awk 'NR > 1 {print $1}' $BASEDIR/samples.txt`; do anvi-profile -i $MAPPING/$sample.bam -c $OUTPUT/contigs.db -o $OUTPUT/$sample --num-threads 18; done
anvi-profile -i $MAPPING/HI.DNA_AS3.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_AS3 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_BS3.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_BS3 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_RS3.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_RS3 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_S1.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_S1 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_S2.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_S2 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_S3.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_S3 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_S4.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_S4 --num-threads 18 --min-contig-length 2500
anvi-profile -i $MAPPING/HI.DNA_S5.bam -c $OUTPUT/contigs.db -o $OUTPUT/HI.DNA_S5 --num-threads 18 --min-contig-length 2500

######
# echo "have things run fine up to here without an error" 
echo "#######"
echo "Message from Graham: Were error messages produced during profiling?" 
echo "#######"
#######

#hmm profiles
echo "Message from Graham: Running HMM" 
anvi-run-hmms -c $OUTPUT/contigs.db --num-threads 18

# annonate genes with NCBI COGs 
#anvi-setup-ncbi-cogs
echo "Message from Graham: Running COGs" 
anvi-run-ncbi-cogs -c $OUTPUT/contigs.db --num-threads 18
# alt option use --sensitive flag to ensure diamond is correctly matching sequences, this is probably unneccessary though

#view contigs stats locally by transferring over contigs.db file
anvi-display-contigs-stats $OUTPUT/contigs.db --report-as-text -o $OUTPUT/contig_stats_anvio5_1.6.txt

#annotate the anvio database CDS (after splitting with http://iubio.bio.indiana.edu/gmod/genogrid/scripts/split_multifasta.pl) for function with InterProScan
anvi-get-aa-sequences-for-gene-calls -c $OUTPUT/contigs.db -o $OUTPUT/protein-sequences.fa 
mkdir -p $OUTPUT/split_sequences/
split_multifasta.pl -i $OUTPUT/protein-sequences.fa -o $OUTPUT/split_sequences/ --seqs_per_file=5000
mkdir -p $OUTPUT/interproscan_output
#run an array of 372 protein multifastas as individual profiles in IPRS and hmmscan for Resfams
#sbatch $OUTPUT/IPRS_gc.sh
#must be run after the annontation
#anvi-import-functions -c $OUTPUT/contigs.db -i $OUTPUT/interpro-output.tsv -p interproscan

#annotate the anvio database contigs for taxonomy with centrifuge
anvi-get-dna-sequences-for-gene-calls -c $OUTPUT/contigs.db -o $OUTPUT/gene-calls.fa
# commented out second line because I don't want to have two separate taxa labels
#centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v $OUTPUT/gene-calls.fa -S $OUTPUT/centrifuge_hits.tsv --threads 18
#anvi-import-taxonomy -c $OUTPUT/contigs.db -i $OUTPUT/centrifuge_report.tsv $OUTPUT/centrifuge_hits.tsv -p centrifuge

#annotate the anvio database contigs for taxonomy with kaiju
$KAIJU/kaiju -t $KAIJUFILES/nodes.dmp -f $KAIJUFILES/kaiju_db_nr.fmi -i $OUTPUT/gene-calls.fa -a greedy -e 5 -m 11 -s 75 -z 18 -o $OUTPUT/gene_calls_nr.out -v
$KAIJU/addTaxonNames -t $KAIJUFILES/nodes.dmp -n $KAIJUFILES/names.dmp -i $OUTPUT/gene_calls_nr.out -o $OUTPUT/gene_calls_nr.names \
              -r superkingdom,phylum,order,class,family,genus,species
$KAIJU/kaijuReport -t $KAIJUFILES/nodes.dmp -n $KAIJUFILES/names.dmp -i $OUTPUT/gene-calls.fa -o $OUTPUT/gene_calls_nr_phylum_Summary.txt -r phylum
$KAIJU/kaijuReport -t $KAIJUFILES/nodes.dmp -n $KAIJUFILES/names.dmp -i $OUTPUT/gene-calls.fa -o $OUTPUT/gene_calls_nr_genus_Summary.txt -r genus
#At this point its not a bad idea to make a copy of your contigs database –just in case.
cp $OUTPUT/contigs.db $OUTPUT/contigs_no_taxa.db
anvi-import-taxonomy-for-genes -i $OUTPUT/gene_calls_nr.names \
                               -c $OUTPUT/contigs.db \
                               -p kaiju --just-do-it

######
# echo "have things run fine up to here without an error" 
echo "#######"
echo "Message from Graham: Did importing kaiju taxa classes appear to work?" 
echo "#######"
#######

#merge everything and bin genomes
anvi-merge $OUTPUT/HI.DNA*/PROFILE.db -o $OUTPUT/SAMPLES-MERGED -c $OUTPUT/contigs.db --sample-name 'LakeHazen2017'
#summarize 
echo "Message from Graham: Summary Beginning, including init-gene-coverage"
anvi-summarize -p $OUTPUT/SAMPLES-MERGED/PROFILE.db -c $OUTPUT/contigs.db \
-o $OUTPUT/SAMPLES-SUMMARY-concoct -C CONCOCT --init-gene-coverages
echo "Message from Graham: Summary Ending"

# need to manually refine bins and resummarize

source deactivate


