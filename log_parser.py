
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import seaborn as sns
import sys
import glob
ORIGINAL_FILE=sys.argv[1]

class BarCodeSeq:
    total=0
    def __init__(self,name,totalreads):
        self.name=name
        self.totalreads=int(totalreads)
        self.sample=None
        self.unique_reads=None
        self.dup_reads=None
        self.unmapped_reads=None
        self.PCR_reads=None

def demultiplexing(title):
    barcode_dict=dict()
    with open('barcodes.bar','r') as barcode_file:
        for line in barcode_file:
            temp=line.split('\t')
            barcode_dict[temp[0]]=temp[1][9:].rstrip()
    with open('log.out','r') as log_file:
        log=log_file.readlines()
    global Barcode_Objects
    Barcode_Objects=[]
    for item in log:
        if item.startswith("Total FastQ records: "):
            temp=item.split(':')
            BarCodeSeq.total=int(temp[1])
            break
    for item in log:
        if item.startswith("FastQ records for barcode"):
            temp=item[26:].split(':')
            Barcode_Objects.append(BarCodeSeq(barcode_dict[temp[0]],temp[1]))


    UNMULTIPLEXED_READS=BarCodeSeq.total-sum([item.totalreads for item in Barcode_Objects])
    fig,ax=plt.subplots(1,1)
    ax.axis('equal')
    ax.set_title('Demultiplexing Output\nfor file {}'.format(title), fontsize=18)
    _,text,__=ax.pie([sum([item.totalreads for item in Barcode_Objects]),UNMULTIPLEXED_READS], explode=(0, 0.3),
       autopct='%1.1f%%', shadow=True, startangle=0,colors=['yellowgreen','lightcoral'],
       labels=['Mapped \nTo\n Barcodes \n({})'.format(sum([item.totalreads for item in Barcode_Objects])),
               'Unmapped \n({})'.format(UNMULTIPLEXED_READS)])
    for item in text:
        item.set_fontsize(14)

    fig.savefig('demultiplexing.pdf')
    mapped_unmapped(title)

def mapped_unmapped(title):
    with open('log.out','r') as log_file:
        for line in log_file:
            if line.startswith('# read'):
                sample=int(line.split(' ')[-1].rstrip())
                mapped=int(next(log_file).split(' ')[-2])
                unmapped=int(next(log_file).split(' ')[-2])
                duplicated=int(next(log_file).split(' ')[-2])
                for item in Barcode_Objects:
                    if item.totalreads==sample:
                        item.unique_reads=mapped
                        item.unmapped_reads=unmapped
                        item.dup_reads=duplicated
                        #assert item.totalreads==item.unique_reads+item.unmapped_reads+item.dup_reads
    ind = np.arange(len(Barcode_Objects))
    Barcode_Objects.sort(key=lambda x:x.totalreads)

    print(ind)
    ind=ind+0.2
    width=0.7
    fig, ax = plt.subplots(figsize=(8,10))

    p1 = ax.barh(ind, [item.unique_reads for item in Barcode_Objects],   width, color='yellowgreen',label="Unique Mapping")
    p2 = ax.barh(ind, [item.dup_reads for item in Barcode_Objects],   width,
              left=np.array([item.unique_reads for item in Barcode_Objects]),
              color='gold',label='Duplicated Mapping')
    p3 = ax.barh(ind, [item.unmapped_reads for item in Barcode_Objects], width, color='lightcoral',label="Unmapped Reads",
             left=np.array([item.unique_reads for item in Barcode_Objects])+np.array([item.dup_reads for item in Barcode_Objects]))
    ax.set_yticks(list(ind+0.4))
    ax.set_yticklabels([item.name for item in Barcode_Objects])

    ax.set_title('Reads Alignment\nfor file {}'.format(title),fontsize=18)
    ax.set_xlabel('Number of Reads Aligned',fontsize=16)
    ax.legend(loc='best')
    fig.savefig('mapping.pdf')



demultiplexing(ORIGINAL_FILE)

def PCR_duplicates(Barcode_Objects=Barcode_Objects):
    log_files=glob.glob('IGV/*.log')
    #global Barcode_Objects
    #print (log_files)
    for log in log_files:
        temp=open(log).readlines()
    #print(temp)
        for item in temp:
        #print(item)
            if item.find('BAR')!=-1:
                bar_line=item
            #print(item)
            else:pass
        barcode=bar_line[bar_line.find('BAR'):bar_line.find('BAR')+5]
        print(barcode)
        for item in temp:
        #print(item)
            if item.find('Unknown')!=-1:
                stats_line=item.split('\t')
                print(stats_line)
            else:pass
        for item in Barcode_Objects:
            if item.name==barcode:
                item.PCR_reads=stats_line[4]

PCR_duplicates()

n_cols=5
n_rows=int(np.ceil(len(Barcode_Objects)/n_cols))
print(int(np.ceil(n_rows)))

fig,ax=plt.subplots(n_cols,n_rows,figsize=(40,40))
max_number_of_reads=[]
for item in Barcode_Objects:
    max_number_of_reads.append(item.totalreads)
    max_number_of_reads.sort(reverse=True)
for n,item in enumerate(ax.flatten()):
    #print(item)
    item.set_title(Barcode_Objects[n].name,fontsize=36)
    item.set_axis_bgcolor('red')
    item.pie([Barcode_Objects[n].unique_reads-int(Barcode_Objects[n].PCR_reads),
        Barcode_Objects[n].dup_reads,
        Barcode_Objects[n].unmapped_reads,
        int(Barcode_Objects[n].PCR_reads)],radius=(Barcode_Objects[n].totalreads)/(max_number_of_reads[0]),colors=['green','blue','red','yellow'])

fig.savefig('final_log.pdf')
