
usage() {
    echo "Usage: analysisWithSpikeIn.sh <inputFile.fq> <runAll or demultiplexOnly>" && exit
}
[[ "$2" != "runAll" && "$2" != "demultiplexOnly" ]] && usage

[ -f $1 ] || { echo Library Fastq file not exist: $1; usage; }

module purge
module load gcc/6.2.0

git clone https://github.com/najoshi/sabre.git
cd sabre
make
cd ..

module purge
module load gcc/6.2.0
module load python/3.6.0
module load fastqc/0.11.3
module load bowtie/1.2.2
module load samtools/0.1.19 ucsc-tools/363 picard/2.8.0 bedtools/2.27.1

#wget https://github.com/BenLangmead/bowtie/archive/v1.1.1.tar.gz >/dev/nul

#tar xvzf v1.1.1.tar.gz >/dev/nul

#export PATH=bowtie-1.1.1:$PATH

#rm -r temp/ IGV/ fastqc/ Jbrowser/ mochiview/ barcodes.bar log.out demultiplexing.pdf 2>/dev/null

sbatch -p short -p priority -n 12 -t 1-0:0:0 --mem 20G -e out.log -o out.log --mail-type=ALL --wrap "python3 analysisWithSpikeIn.py ${1} ${2}"
 
#python3 analysis.py ${1}
