
module load seq/fastx/0.0.13
module load dev/python/3.4.2
module load seq/bowtie/1.1.1
module load seq/samtools/1.2

bsub -q priority -n 4 -W 24:00 -o log.out  "python3.4 analysis.py ${1}"
