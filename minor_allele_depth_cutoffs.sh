#bin/bash

#HEADER FOR SUBMITTED SCRIPTS
#SBATCH --job-name=Cutoffs
#SBATCH --partition=batch 
#SBATCH --ntasks=1
#SBATCH  --nodes=1 
#SBATCH  --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --export=NONE
#SBATCH --mem=20gb
#SBATCH --mail-user=drt83172@uga.edus
#SBATCH --mail-type=END,FAIL
#SBATCH --output=%x_%j.out 
#SBATCH --error=%x_%j.err 



# cd $PBS_O_WORKDIR #to use the directory from which the job is submitted as the working directory (where to find input files or binaries)

## following lines are optional to record info in stdout file
# echo
# echo "Job ID: $PBS_JOBID"
# echo "Queue:  $PBS_QUEUE"
# echo "Cores:  $PBS_NP"
# echo "Nodes:  $(cat $PBS_NODEFILE | sort -u | tr '\n' ' ')"
# echo "mpirun: $(which mpirun)"
# echo

parents3000_all = "/scratch/drt83172/Wallace_lab/TallFescue/Data/VCFs/parents3000_all.vcf"
UGA_149001_FlexSeqResults = "/scratch/drt83172/Wallace_lab/TallFescue/Data/VCFs/UGA_149001_FlexSeqResults.vcf"

# Use this to create a list of values to import to R based on one data measurment. Allows us to create cutoffs for depth and MAC. 
# awk 'BEGIN{found=0;}/DP=/{if(!found){found=1;$0=substr($0,index($0, "DP=")+3);}}/\;/{if(found){found=2;$0=substr($0,0,index($0,";")-1);}}{ if(found){print;if(found==2)found=0;}}'

# The below command will help find depth cut offs much better than the above one
# --vcf UGA_149001_FlexSeqResults.vcf --mac 20  --site-mean-depth

ml BCFtools/1.10.2-GCC-8.3.0
module load VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0


# Past filters for minor allel frequency was to include sites greater than or equal to .0441 for parental data
#Code on minor allel cutoff and depth cut off will go here
# Keeps 1305
vcftools --vcf UGA_149001_FlexSeqResults.vcf --min-meanDP 20 --max-meanDP 600 --mac 20 --out UGA_149001_FlexSeqResults_filtered

# Keeps 1213, this is what I used 
vcftools --gzvcf UGA_149001_FlexSeqResults.vcf.gz --min-meanDP 45 --max-meanDP 600 --mac 20 --max-missing .25 --recode   --recode-INFO-all --out UGA_149001_FlexSeqResults_filtered

# This command will find the snps we are using from the flex seq data in the parental data
# First make an index
bgzip UGA_149001_FlexSeqResults_filtered.recode.vcf
bcftools index UGA_149001_FlexSeqResults_filtered.recode.vcf.gz
gunzip UGA_149001_FlexSeqResults_filtered.recode.vcf.gz
bcftools index parents.bcf
bcftools view --regions-file UGA_149001_FlexSeqResults_filtered.recode.vcf --output-type b --output-file parents1500.bcf parents.bcf


# The above command did not retrive all the needed snp positions. I will now make a list of missing psitions and extract them.
bcftools sort --output-type b --output-file parents1500_sorted.bcf parents1500.bcf
bcftools sort --output-type v --output-file UGA_149001_FlexSeqResults_filtered.recode_sorted.vcf UGA_149001_FlexSeqResults_filtered.recode.vcf
#Gives list of all chromosomes and positions for flex data
sed -n '1864,3077p' UGA_149001_FlexSeqResults_filtered.recode_sorted.vcf | cut -f 1,2 > list_all.txt
bcftools view parents1500_sorted.bcf | sed -n '2209119,2210317p' | cut -f 1,2 > list_have.txt
diff list_all.txt list_have.txt | grep "<" | sed s'/<//'g | sed s'/ //'g> list_need.txt


bcftools view --targets-file list_need.txt --output-type b --output-file parents1500_missing.bcf parents.bcf
bcftools index parents1500_missing.bcf
bcftools index parents1500.bcf
bcftools concat -a --output-type v --output parents1500_all.vcf parents1500_missing.bcf parents1500.bcf
bcftools sort --output-type v --output-file parents1500_all_sorted.vcf parents1500_all.bcf

run_pipeline.pl -SortGenotypeFilePlugin -inputFile ~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/parents1500_all.vcf -outputFile ~/Documents/Tall_fescue/Genotype_Data/2021_01_22_FescueFlexSeqGenos/parents1500_all_sorted.vcf
