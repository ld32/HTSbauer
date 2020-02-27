#author luis soares

import sys
import os
import glob
import subprocess
from distutils.dir_util import copy_tree

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]


CURRENT_PATH=os.getcwd()
ORIGINAL_FILE=sys.argv[1]
action=sys.argv[2]

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
			barcodes.append(line.rstrip().split())
	return (demult_par,bowtie_par,macs2_par,barcodes)

def create_barcode_files(barcodes):
	needed=[]
	for item in barcodes:
		needed.append(item[1].upper())
		needed.append(item[2].upper())
	needed=set(needed)
	bar=open('barcodes.bar','w')
	for item in needed:
		bar.write(ALL_BARCODES[item]+'\t'+'temp_bar_'+item+'\n')
	bar.close()

def demultiplexing(ORIGINAL_FILE):
	exit_code=subprocess.call("sabre/sabre se -m 1 -f {0} -b barcodes.bar -u temp_bar_unmatched".format(ORIGINAL_FILE), shell=True)  
	if exit_code!=0:
		sys.exit()
	files=os.listdir()
	for item in files:
		if 'temp_bar_' in item:
			os.rename(item,CURRENT_PATH+'/temp/'+item + ".fq")

def fastqc():
	if os.path.isfile('qcIsDone.txt'):
		return
	files=glob.glob('temp/*.fq')
	files.remove('temp/temp_bar_unmatched.fq')
	temp=chunks(files,10)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				processes.append(subprocess.Popen("fastqc {0}".format(item),shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:
		pass
	os.makedirs('fastqc',exist_ok=True)
	files=glob.glob('temp/*.html')
	files.extend(glob.glob('temp/*.zip'))
	for item in files:
		os.rename(item,CURRENT_PATH+'/fastqc/'+os.path.basename(item))
	os.system("touch qcIsDone.txt")

def align_spombe(): 
	files=glob.glob('temp/*.fq')
	files.remove('temp/temp_bar_unmatched.fq')
	temp=chunks(files,1)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				print("working on " + item)
				#if(os.stat(item).st_size != 0):
				processes.append(subprocess.Popen("bowtie -S -p 10 indexes/Spombe {0}.fq | samtools import indexes/Spombe.fa.fai - {0}.bam".format(item[:-3]),shell=True))
			exit_codes=[p.wait() for p in processes]
			
	except StopIteration:pass

def split_bams(): 
	files=glob.glob('temp/*.fq')
	files.remove('temp/temp_bar_unmatched.fq')
	temp=chunks(files,5)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				if(os.stat(item).st_size != 0):
					processes.append(subprocess.Popen("samtools view -b -F4 {0}.bam -o {0}_align.bam".format(item[:-3]),shell=True))
					processes.append(subprocess.Popen("samtools view -b -f4 {0}.bam -o {0}_unalign.bam".format(item[:-3]),shell=True))
			exit_codes=[p.wait() for p in processes]
			
	except StopIteration:pass

def convert_bam_to_fq(): 
	files=glob.glob('temp/*.fq')
	files.remove('temp/temp_bar_unmatched.fq')
	temp=chunks(files,5)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				if(os.stat(item).st_size != 0):
					processes.append(subprocess.Popen("bedtools bamtofastq -i temp/{0}_align.bam -fq temp1/{0}_aligned.fq".format(os.path.basename(item)[:-3]),shell=True))
					processes.append(subprocess.Popen("bedtools bamtofastq -i temp/{0}_unalign.bam -fq temp1/{0}_unalign.fq".format(os.path.basename(item)[:-3]),shell=True))
			exit_codes=[p.wait() for p in processes]		
	except StopIteration:pass

def align_genome_aligned():
	os.system("rm alignedBowtie.log 2>/dev/null")
	files=glob.glob('temp1/*aligned.fq')
	temp=chunks(files,1)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				if(os.stat(item).st_size != 0):
					processes.append(subprocess.Popen("echo working on:{0} >>alignedBowtie.log; bowtie -S -p 10 indexes/genome {0}_aligned.fq 2>>alignedBowtie.log | samtools view -bS -F 4 - > {0}pom_unalign.bam".format(item[:-11]),shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass

def align_genome_unaligned():
	os.system("rm unalignBowtie.log 2>/dev/null") 
	files=glob.glob('temp1/*unalign.fq')
	temp=chunks(files,1)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				if(os.stat(item).st_size != 0):
					processes.append(subprocess.Popen("echo working on:{0} >>unalignBowtie.log; bowtie -S -p 10 -m 1 indexes/genome temp1/{0}_unalign.fq 2>>unalignBowtie.log | samtools view -bS -F 4 - > temp2/{0}.bam".format(os.path.basename(item)[:-11]),shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass


def sam_tools():
	files=glob.glob('temp2/*.bam')
	temp=chunks(files,10)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				processes.append(subprocess.Popen("samtools sort {0}.bam {0}_sorted".format(item[:-4]),shell=True))
			exit_codes=[p.wait() for p in processes]
			processes=[]
			for item in proc_files:
				processes.append(subprocess.Popen("samtools index {0}_sorted.bam".format(item[:-4]),shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass
	os.makedirs('IGV',exist_ok=True)
	files=glob.glob('temp2/*sorted*')
	for item in files:
		os.rename(item,CURRENT_PATH+'/IGV/'+os.path.basename(item))
		
def macs2():
	# This is to only run one copy of macs2 to avoid race condition
	if not os.path.exists(os.environ['HOME'] + "/.python-eggs/MACS2-2.1.1.20160309-py2.7-linux-x86_64.egg-tmp/MACS2"):
		print("python egg cache folder not exist")
		temp=chunks(BARCODES,1)
		try:
			while True:
				processes=[]
				proc_files=next(temp)	
				for item in proc_files:
					name=item[0]
					input_bar=item[1]
					sample_bar=item[2]
					processes.append(subprocess.Popen("module load python/2.7.12 macs2/2.1.1.20160309 2>/dev/null; macs2 callpeak -t temp2/temp_bar_"+sample_bar+".bam -c temp2/temp_bar_"+input_bar+".bam -f BAM -g 12100000 -n temp2/"+name+" -B -q 0.01 --nomodel --extsize 150 --SPMR",shell=True))
				exit_codes=[p.wait() for p in processes]
				break
		except StopIteration:pass
		
	# After the python egg cache built ready, macs2 can run in parallel
	temp=chunks(BARCODES,10)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				name=item[0]
				input_bar=item[1]
				sample_bar=item[2]
				processes.append(subprocess.Popen("module load python/2.7.12 macs2/2.1.1.20160309 2>/dev/null; macs2 callpeak -t temp2/temp_bar_"+sample_bar+".bam -c temp2/temp_bar_"+input_bar+".bam -f BAM -g 12100000 -n temp2/"+name+" -B -q 0.01 --nomodel --extsize 150 --SPMR",shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass

def wig():
	#subprocess.call('module purge',shell=True)
	#subprocess.call('module load dev/python/3.4.2',shell=True)
	files=glob.glob('temp2/*')
	input_files=[]
	for item in files:
		if "treat_pileup" in item:
			input_files.append(item)
	temp=chunks(input_files,10)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				processes.append(subprocess.Popen('python3 BDGtoWIG.py -i '+item,shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass

def batch_wig():
	os.makedirs('mochiview',exist_ok=True)
	files=glob.glob('temp2/*.wig')
	for item in files:
		os.rename(item,CURRENT_PATH+'/mochiview/'+os.path.basename(item))
		
def bigwig():
	files=glob.glob('mochiview/*.wig')
	#subprocess.call('module load seq/UCSC-tools',shell=True)
	temp=chunks(files,11)
	try:
		while True:
			processes=[]
			proc_files=next(temp)	
			for item in proc_files:
				processes.append(subprocess.Popen("wigToBigWig {0} sizes {0}.bwig".format(item),shell=True))
			exit_codes=[p.wait() for p in processes]
	except StopIteration:pass

def batch_bwig():
	os.makedirs('Jbrowser',exist_ok=True)
	files=glob.glob('mochiview/*.bwig')
	for item in files:
		os.rename(item,CURRENT_PATH+'/Jbrowser/'+os.path.basename(item))
		
def duplicate_remove():
	files=glob.glob('IGV/*.bam')
	for item in files:
		name=item[9:18]
		subprocess.call("java -Xmx1G -jar $PICARD/picard-2.8.0.jar MarkDuplicates INPUT={0} OUTPUT={0}nodup.bam METRICS_FILE=IGV/{1}.log REMOVE_DUPLICATES=true ASSUME_SORTED=false".format(item,name),shell=True)

DEMULT_PAR,BOWTIE_PAR,MACS2_PAR,BARCODES=parse_setup('setup.cfg')
print(CURRENT_PATH)
print(DEMULT_PAR,BOWTIE_PAR,MACS2_PAR,BARCODES)
create_barcode_files(BARCODES)
os.makedirs('temp', exist_ok=True)
os.makedirs('temp1', exist_ok=True)
os.makedirs('temp2', exist_ok=True)
demultiplexing(ORIGINAL_FILE)

if action == "demultiplexOnly": 
    quit()

fastqc()
align_spombe()
split_bams()
convert_bam_to_fq()
align_genome_aligned()
align_genome_unaligned()
sam_tools()
macs2()
wig()
batch_wig()
bigwig()
batch_bwig()
duplicate_remove()
 
subprocess.call('module load python/2.7.12 2>/dev/null; python log_parserWithSpikeIn.py {}'.format(ORIGINAL_FILE),shell=True)
