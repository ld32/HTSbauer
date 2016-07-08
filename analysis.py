#author luis soares

import sys
import os
import glob
import subprocess
from distutils.dir_util import copy_tree

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


CURRENT_PATH=os.getcwd()
ORIGINAL_FILE=sys.argv[1]
ALL_BARCODES={'BAR13':'ATATAGGA',
			 'BAR14':'AACCGTGT',
			'BAR15':'AGGTCAGT',
			'BAR16':'CTCTGTCT',
			'BAR17':'CCATACAC',
			'BAR18':'CGCATTAA',
			'BAR19':'GTCTACAT',
			'BAR20':'GAGTTAAC',
			'BAR21':'GCAGCCTC',
			'BAR22':'TCGCGTAC',
			'BAR23':'TATACCGT',
			'BAR24':'TGCGGTTA',
			'BAR28':'TTGAGTGT',
			'BAR32':'TGACGCAT',
			'BAR36':'TGATCCGA',
			'BAR40':'GTGGGATA',
			'BAR44':'TATCTCCG',
			'BAR48':'TGAGAGTG',
			'BAR52':'TTCTGATG',
			'BAR56':'TCATTAGG',
			'BAR60':'AGAACACC',
			'BAR25':'AACACCTAC',
			'BAR26':'CCTTTACAG',
			'BAR27':'GGTCCTTGA',
			'BAR29':'ACTAACTGC',
			'BAR30':'CAGGAGGCG',
			'BAR31':'GTTGTCCCA',
			'BAR33':'ATCGCCAGC',
			'BAR34':'CATTCCAAG',
			'BAR35':'GCAAGTAGA',
			'BAR37':'ACGTAGCTC',
			'BAR38':'CGAACTGTG',
			'BAR39':'TAGCTAGTA',
			'BAR41':'ATCCTATTC',
			'BAR42':'CGGACGTGG',
			'BAR43':'GCGTTTCGA',
			'BAR45':'ACAGTGCAC',
			'BAR46':'CACAGTTGG',
			'BAR47':'GTGACTACA',
			'BAR49':'AATGCTGAC',
			'BAR50':'CCGTCTGAG',
			'BAR51':'GGCAGACGA',
			'BAR53':'AGTAGTGGC',
			'BAR54':'CTAGTCATG',
			'BAR55':'GACACTCTA',
			'BAR57':'TCCAGCCTC',
			'BAR58':'CTAGATTCG',
			'BAR59':'GAACGCTGA'}

			

def parse_setup(setup):
	with open(setup) as init:
		demult_par=init.readline().split(':')[1].rstrip()
		bowtie_par=init.readline().split(':')[1].rstrip()
		macs2_par=init.readline().split(':')[1].rstrip()
		barcodes=[]
		for line in init:
			barcodes.append(line.rstrip().split('\t'))
	return (demult_par,bowtie_par,macs2_par,barcodes)
	
def create_barcode_files(barcodes):
	needed=[]
	for item in barcodes:
		needed.append(item[1].upper())
		needed.append(item[2].upper())
	needed=set(needed)
	size8=open('size8.bar','w')
	size9=open('size9.bar','w')
	for item in needed:
		if len(ALL_BARCODES[item])==8:
			size8.write(item+'\t'+ALL_BARCODES[item]+'\n')
		elif len(ALL_BARCODES[item])==9:
			size9.write(item+'\t'+ALL_BARCODES[item]+'\n')
	size8.close()
	size9.close()
	
def demultiplexing(ORIGINAL_FILE):
	exit_code=subprocess.call("zcat {0} | fastx_barcode_splitter.pl --bcfile size8.bar  --prefix temp_bar_ --bol --mismatches {1}".format(ORIGINAL_FILE,DEMULT_PAR), shell=True)  
	if exit_code!=0:
		sys.exit()
	os.rename('temp_bar_unmatched','unmatched')
	exit_code=subprocess.call("cat {0} | fastx_barcode_splitter.pl --bcfile size9.bar  --prefix temp_bar_ --bol --mismatches {1}".format('unmatched',DEMULT_PAR), shell=True)
	if exit_code!=0:
		sys.exit()
	os.remove('unmatched')
	files=os.listdir()
	for item in files:
		if 'temp_bar_' in item:
			os.rename(item,CURRENT_PATH+'/temp/'+item)
			
def trim_fastq():
	files=glob.glob('temp/*')
	files.remove('temp/temp_bar_unmatched')
	for item in files:
		subprocess.call("python3.4 fastq_trimmer.py -f {0} -i {1} -o {2}_trim".format(str(len(ALL_BARCODES[item[-5:]])+1),item,item),shell=True)
		os.remove(item)

def fastqc():
	files=glob.glob('temp/*')
	files.remove('temp/temp_bar_unmatched')
	for item in files:
		subprocess.call("fastq {0}".format(item),shell=True)
		
def batch_fastqc():
	os.makedirs('fastqc',exist_ok=True)
	files=glob.glob('temp/*.html')
	files.extend(glob.glob('temp/*.zip)'))
	for item in files:
		os.rename(item,CURRENT_PATH+'/mochiview/'+os.path.basename(item))

def bowtie_align():
	files=glob.glob('temp/*')
	files.remove('temp/temp_bar_unmatched')
	os.makedirs('indexes')
	copy_tree('/n/data2/hms/bcmp/buratowski/indexes', CURRENT_PATH+'/indexes')
	for item in files:
		subprocess.call("bowtie -S -p 8 -m 1 genome {0} {0}.sam".format(item),shell=True)
		os.remove(item)
		
def sam_tools():
	files=glob.glob('temp/*')
	files.remove('temp/temp_bar_unmatched')
	for item in files:
		subprocess.call("samtools import genome.fa.fai {0} {1}bam".format(item,item[:-3]),shell=True)
		os.remove(item)
		subprocess.call("samtools sort {0}bam {0}sorted".format(item[:-3]),shell=True)
		os.remove(item[:-3]+'bam')
		subprocess.call("samtools index {0}sorted.bam".format(item[:-3]),shell=True)
		
def macs2():
	subprocess.call('module purge',shell=True)
	subprocess.call('module load seq/macs/2.1.0',shell=True)
	subprocess.call('module load dev/python/2.7.6',shell=True)
	for item in BARCODES:
		name=item[0]
		input_bar=item[1]
		sample_bar=item[2]
		subprocess.call("macs2 callpeak -t temp/temp_bar_"+sample_bar+"_trim.sorted.bam -c temp/temp_bar_"+input_bar+"_trim.sorted.bam -f BAM -g 12100000 -n temp/"+name+" -B -q 0.01 --nomodel --extsize 150 --SPMR",shell=True)
		
def wig():
	subprocess.call('module purge',shell=True)
	subprocess.call('module load dev/python/3.4.2',shell=True)
	files=glob.glob('temp/*')
	input_files=[]
	for item in files:
		if "treat_pileup" in item:
			input_files.append(item)
	for item in input_files:
		subprocess.call('python3.4 BDGtoWIG.py -i '+item,shell=True)

def batch_wig():
	os.makedirs('mochiview',exist_ok=True)
	files=glob.glob('temp/*.wig')
	for item in files:
		os.rename(item,CURRENT_PATH+'/mochiview/'+os.path.basename(item))

DEMULT_PAR,BOWTIE_PAR,MACS2_PAR,BARCODES=parse_setup('setup.cfg')
#print(CURRENT_PATH)
#print(DEMULT_PAR,BOWTIE_PAR,MACS2_PAR,BARCODES)
create_barcode_files(BARCODES)
os.makedirs('temp', exist_ok=True)
demultiplexing(ORIGINAL_FILE)
trim_fastq()
fastqc()
batch_fastqc()
bowtie_align()
sam_tools()
macs2()
wig()
batch_wig()

subprocess.call('python2.7 log_parser.py {}'.format(ORIGINAL_FILE),shell=True)
