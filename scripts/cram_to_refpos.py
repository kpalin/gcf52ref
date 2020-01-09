# wdecoster

import pysam
import numpy as np
from typing import Tuple,List,Dict
from argparse import ArgumentParser
import sys



from functools import lru_cache
class Reference(object):
    def __init__(self,fastafile):
        self.fastafile=fastafile
        import pysam
        self._fasta = pysam.FastaFile(self.fastafile)

    @lru_cache(maxsize=10000)
    def fetch(self,chrom,start,end):
        return self._fasta.fetch(chrom,start,end)


def main():
    import sys
    args = get_args()
    if not args.no_header:
        print("\t".join(['#chromosome', 'strand', 'start',
                         'end', 'quality', 'read_name', 'modification',"sequence"]))

    if args.reference_fasta is not None:
        fasta = Reference(args.reference_fasta)
    else:
        class DummyReference(object):
            def fetch(self,*args,**kwargs):
                return "N"
        fasta = DummyReference()

    cram = pysam.AlignmentFile(args.cram, "rc")
    outstrm = sys.stdout
    for read in cram.fetch(reference=args.chrom, start=args.start, end=args.end):
        if read.has_tag('MM') and not read.is_supplementary and not read.is_secondary:
            mod, positions, qualities = get_modified_reference_positions(read)
            for pos, qual in zip(positions, qualities):
                if pos is not None:
                    outstrm.write("{chrom}\t{strand}\t{st}\t{en}\t{qual}\t{rname}\t{mod}\t{ref}\n".format(
                        chrom=read.reference_name,
                        strand='-' if read.is_reverse else '+',
                        st=pos,
                        en=pos,
                        qual=qual,
                        rname=read.query_name,
                        mod=mod,
                        ref=fasta.fetch(read.reference_name,pos-3,pos+4)))

def get_likelihoods(read) -> Dict[str,np.ndarray]:
    "If the read has ml:Z: tag, return the modification likelihoods as numpy array for each base of the read"
    import numpy as np
    
    likelihoods = dict()
    mls = read.get_tag("ml").rstrip(";").split(";")
    for ml in mls:
        basemod,ml=ml.split(",",1)
        ml = np.frombuffer(ml.encode("ascii"),dtype="|S2")
        ml = np.fromiter((int(x,16) for x in ml),dtype="u1") 
        likelihoods[basemod] = ml
        assert len(ml) == read.infer_read_length()

    return likelihoods


def get_modified_reference_positions(read) -> Tuple[str,np.ndarray,List[int]]:
    import numpy as np
    basemod = read.get_tag('MM').split(',', 1)[0]
    if '-' in basemod:
        sys.exit("ERROR: modifications on negative strand currently unsupported.")
    base, mod = basemod.split('+')
    deltas = [int(i) for i in read.get_tag('MM').split(',')[1:]]
    quals = [ord(i) - 33 for i in read.get_tag('MP')]
    locations = np.cumsum(np.array(deltas,dtype=int)+1) - 1
 
    read_sequence = np.frombuffer(read.get_forward_sequence().encode("ascii"),dtype="|S1")
    base_index, = np.where(read_sequence == base.encode("ascii") )

    modified_bases = base_index[locations]

    refpos = read.get_reference_positions(full_length=True)
    #likelihoods = get_likelihoods(read)[basemod][modified_bases]
    if read.is_reverse:
        return basemod, [refpos[-i-1] for i in modified_bases], quals[::-1]
    else:
        return basemod, [refpos[i] for i in modified_bases], quals


def get_args():
    parser = ArgumentParser(description="convert ont-cram file with MM/MP tags to tsv")
    parser.add_argument("cram", help="aligned cram file containing MM/MP tags")
    parser.add_argument("--chrom", help="chromosome to limit calls to.", default=None)
    parser.add_argument("--start", help="start coordinates to limit calls from",type=int)
    parser.add_argument("--end", help="end coordinates to limit calls from",type=int)
    parser.add_argument("--reference_fasta", help="indexed fasta file for reporting sequence context.", default=None)
    parser.add_argument("--no_header",help="Don't output header line",action="store_true",default=False)
    return parser.parse_args()


if __name__ == '__main__':
    main()
