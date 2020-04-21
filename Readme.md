This is forked from: LuisSoares/HTSbauer

Added function to take care of spikein library. 

# Instructions to use HTS pipeline

This pipeline works for raw files coming from the Bauer center sequencing facility.

sh analysisWithSpikeIn.sh HJ2_NoIndex.R1.fastq.gz runAll

or

sh analysisWithSpikeIn.sh HJ2_NoIndex.R1.fastq.gz demultiplexOnly


Note: this workflow does not normallize the data, only calcuate the spombe/genome ratios.


For future reference, here is manual version of the workflow. 

After demultiplexing and trimming the fastq files, you need to align each barcode file to pombe, the indexes in the lab folder now include the pombe genome under the name Spombe, so you need to run:
 
bowtie –S –p 8 Spombe <fq file> <sam_file>.sam    #notice that you don’t really need the –m 1
 
from the output you should take the total number of reads, the number of aligned reads and the number of unaligned reads (this numbers are not need anywhere but they are a good sanity check for the later steps), if you don’t want the output in the mail just pass the –o <log_file_name> into the bsub command.
 
Then convert the sam file to bam again using Spombe instead of genome:
 
samtools import Spombe.fa.fai <sam_file> <bam_file>.bam
 
Next you need to separate the bam file into aligned and unaligned reads (for the pombe genome)
 
samtools view –h –b -F4 <bam_file> > <bam_file_aligned>  # the flag –F4 gives you just the aligned alignments
samtools view –h –b –f4 <bam_file> > <bam_file_unaligned>  # the flag –f4 gives you just the unaligned alignments
 
you then need to import BEDTOOLS, and run the following for both the aligned and unaligned bam files from the previous step
 
bedtools bamtofastq –i <bam_file_aligned> –fq <bam_file_aligned>.fq
bedtools bamtofastq –i <bam_file_unaligned> –fq <bam_file_unaligned>.fq
 
So now you have two need fastq files, which you are going to use for bowtie with cerevisiae this time:
 
bowtie –S –p 8 genome <bam_file_aligned.fq file> <sam_file>.sam    #notice that you don’t really need the –m 1
bowtie –S –p 8 –m 1 genome <bam_file_unaligned.fq file> <sam_file>.sam    #notice that you need the –m 1
 
from the output of the first bowtie you will get the aligned number (those will be the ambiguous reads) and the non-aligned reads (those are the
unambiguously pombe reads)
from the output of the second bowtie the sum of aligned and multialigned are the unambiguously  cerevisiae reads, the unaligned are the contamination reads (the sam file output of the second bowtie you continue with the normal pipeline, sort, macs, and wig, that will be the non-normalized cerevisiae tracks).
Now sum the unambiguously pombe and unambiguously cerevisiae reads, those are your total reads, then for the inputs and IPs, calculate pombe/total,
then for just the IPs divide the Input ratio by the IP ratio, this will be the normalization factor than can be used for each track (this is how they do it in the Winston lab and more or less how it is described in the paper: Biological chromodynamics: a general method for measuring protein occupancy across the genome by calibrating Chip-Seq.), I am actually using the square root of the normalization factor as a normalization factor since it looks much better (I will try to come up with a possible explanation for this, maybe even the square root is not appropriate and maybe it should be adjusted by some statistical function that describes the probability of a read be sequenced). 
