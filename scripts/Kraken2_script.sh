# Prefetch the SRA files from NCBI
for line in $(cat sra.txt); do echo "Prefetching $line"; prefetch $line; echo "$line is done"; done

# One didn't work
prefetch SRR8146944

# Now download the fastq files
for file in ./*/*.sra; do echo "Getting fastq $file"; fasterq-dump $file --split-files --threads 8 --progress --details --temp ./tmp/ --outdir ./out/ ; echo "fastq $file is downloaded"; done

# Now download database
ml kraken2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz 
tar -xvf k2_core_nt_20251015.tar.gz

# Checking quality with FastQC
ml fastqc
for file in ./out/*.fastq;
do
fastqc -o ./fastqc_out/ $file;
done

# In order to trim with cutadapt, install the python wheel
avail_wheels cutadapt
module avail python
ml python
virtualenv --no-download ./
source ./bin/activate
pip install --no-index --upgrade pip
pip install cutadapt --no-index

# Now trim
mkdir cutadapt_out
for file in ./out/*_1.fastq; do
base=$(basename "$file" "_1.fastq")
cutadapt -q 20,20 -o ./cutadapt_out/${base}_1_trimmed.fastq -p ./cutadapt_out/${base}_2_trimmed.fastq ./out/${base}_1.fastq ./out/${base}_2.fastq;
done

deactivate

# Check again with FastQC
ml fastqc
for file in ./*.fastq;
do
fastqc -o ./fastqc_trimmed_out/ $file;
done

# Script for classifying reads with kraken2
#!/bin/bash
#SBATCH --job-name=classify        # Job name
#SBATCH --nodes=1                 # Number of nodes
#SBATCH --cpus-per-task=16    # CPUs per task
#SBATCH --mem=64G                  # Memory per node 
#SBATCH --time=12:00:00           # Time limit
#SBATCH --output=output_%j.txt   # Output file
#SBATCH --error=error_%j.txt     # Error file

ml kraken2

for file in ./cutadapt_out/*_1_trimmed.fastq;
do
base=$(basename "$file" "_1_trimmed.fastq")
kraken2 --db ./db/ --confidence 0.15 --output "${base}.kraken" --report "${base}.report" --paired --use-names "./cutadapt_out/${base}_1_trimmed.fastq" "./cutadapt_out/${base}_2_trimmed.fastq";
done

# End script #

# Abundance estimation with bracken

ml bracken

for file in *.report;
do
base=$(basename "$file" ".report")
bracken -d ../db/ -i $file -o ${base}.bracken
echo "Abundance of $base estimated";
done

# Install python environment of kraken-biom

ml python
virtualenv --no-download kraken-biom
source kraken-biom/bin/activate
pip install --no-index --upgrade pip
pip install kraken-biom

# Now converting to BIOM format 
kraken-biom SRR8146935_bracken_species.report SRR8146938_bracken_species.report SRR8146951_bracken_species.report SRR8146954_bracken_species.report SRR8146936_bracken_species.report SRR8146944_bracken_species.report SRR8146952_bracken_species.report SRR8146956_bracken_species.report --fmt json

# End




