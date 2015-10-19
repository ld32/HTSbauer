__author__ = 'Luis'
import numpy as np
from collections import OrderedDict
import argparse

parser = argparse.ArgumentParser(description='Split FASTQ file according to barcodes')
parser.add_argument('-i', dest='input_file', help='input file',nargs=1,required=True)
args = parser.parse_args()

input_file=args.input_file[0]

chromosome_dict=OrderedDict()
list_of_chromosome=['I',
                    'II',
                    'III',
                    'IV',
                    'V',
                    'VI',
                    'VII',
                    'VIII',
                    'IX',
                    'X',
                    'XI',
                    'XII',
                    'XIII',
                    'XIV',
                    'XV',
                    'XVI']

size_of_chromosomes=[230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066]

for n in range(len(list_of_chromosome)):
    chromosome_dict[list_of_chromosome[n]]=np.zeros(size_of_chromosomes[n])


with open(input_file) as data:
    for line in data:
        if 'MT' in line:continue
        line_list=line.rstrip().split('\t')
        chromosome_dict[line_list[0]][int(line_list[1]):int(line_list[2])]=float(line_list[3])


with open(input_file+'.wig','w') as output:
    output.write('track\n')
    for chrm in chromosome_dict.keys():
            output.write('variableStep chrom=chr{}\n'.format(chrm))
            n=0
            while n<len(chromosome_dict[chrm]):
                output.write(str(n+1)+'\t'+str(chromosome_dict[chrm][n])+'\n')
                n+=1
