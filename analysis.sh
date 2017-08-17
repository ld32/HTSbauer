

module load seq/macs/2.1.0
module load seq/fastqc/0.11.3
module load seq/fastx/0.0.13
module load dev/python/3.4.2
module load seq/bowtie/1.1.1
module load seq/samtools/1.2
module load dev/python/2.7.6
module load dev/gcc-5.2.0
module load seq/picard/1.138
module load seq/UCSC-tools

git clone https://github.com/najoshi/sabre.git
cd sabre
make
cd ..

bsub -q priority -n 12 -W 24:00 -N -oo log.out  "python3.4 analysis.py ${1}"


