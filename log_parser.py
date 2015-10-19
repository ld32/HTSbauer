__author__ = 'Luis'

class BarCodeSeq:
    def __init__(self,name,totalreads):
        self.name=name
        self.totalreads=totalreads
        self.sample=None
        self.unique_reads=None
        self.dup_reads=None
        self.unmapped_reads=None


log=open('log.out','r')
