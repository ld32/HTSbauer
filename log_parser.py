__author__ = 'Luis'
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class BarCodeSeq:
    total=0
    def __init__(self,name,totalreads):
        self.name=name
        self.totalreads=int(totalreads)
        self.sample=None
        self.unique_reads=None
        self.dup_reads=None
        self.unmapped_reads=None

def demultiplexing(title):
    with open('log.out','r') as log_file:
        log=log_file.readlines()
    global Barcode_Objects
    Barcode_Objects=[]
    for item in log:
        if item.startswith("BAR"):
            temp=item.split('\t')
            Barcode_Objects.append(BarCodeSeq(temp[0],temp[1]))

    for item in log:
        if item.startswith("total"):
            temp=item.split('\t')
            BarCodeSeq.total=int(temp[1])
            break

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

    plt.savefig('demultiplexing.pdf')
    plt.close()
    mapped_unmapped(title)

def mapped_unmapped(title):
    with open('log.out','r') as log_file:
        for line in log_file:
            if line.startswith('# read'):
                sample=int(line.split(' ')[-1].rstrip())
                mapped=int(log_file.readline().split(' ')[-2])
                unmapped=int(log_file.readline().split(' ')[-2])
                duplicated=int(log_file.readline().split(' ')[-2])
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

    p1 = ax.barh(ind, [item.unique_reads for item in Barcode_Objects],   width, color='yellowgreen',label="Unique Reads")
    p2 = ax.barh(ind, [item.dup_reads for item in Barcode_Objects],   width,
              left=np.array([item.unique_reads for item in Barcode_Objects]),
              color='gold',label='Duplicated Reads')
    p3 = ax.barh(ind, [item.unmapped_reads for item in Barcode_Objects], width, color='lightcoral',label="Unmapped Reads",
             left=np.array([item.unique_reads for item in Barcode_Objects])+np.array([item.dup_reads for item in Barcode_Objects]))
    plt.yticks(ind+0.4,[item.name for item in Barcode_Objects],rotation='horizontal')
    plt.title('Reads Alignment\nfor file {}'.format(title),fontsize=18)
    plt.xlabel('Number of Reads Aligned',fontsize=16)
    plt.legend(loc='best')
    plt.savefig('mapping.pdf')

