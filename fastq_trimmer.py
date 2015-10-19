#trimms fastq file sequences, removing barcode size
#python3.3

import argparse

parser = argparse.ArgumentParser(description='Split FASTQ file according to barcodes')
parser.add_argument('-f', dest='first', help='first base to keep',nargs=1,required=True)
parser.add_argument('-i', dest='input_file', help='input file',nargs=1,required=True)
parser.add_argument('-o', dest='output_file', help='output file prefix',nargs=1,required=True)
args = parser.parse_args()

first=int(args.first[0])
input_file=args.input_file[0]
output_file=args.output_file[0]


def create_output_files(input_file,output_file,first):
    starting_file=open(input_file,'r')
    out_file=open(output_file,'w')
    header='start'
    while header:
        header=starting_file.readline()
        sequence=starting_file.readline()
        header=header+sequence[first:]
        header=header+starting_file.readline()
        header=header+starting_file.readline()[first:]
        out_file.write(header)
    out_file.close()

if __name__=='__main__':
    create_output_files(input_file,output_file,first)