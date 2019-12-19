# wdecoster

import pysam
import numpy as np
from argparse import ArgumentParser
import sys


def main():
    args = get_args()
    print("\t".join(['chromosome', 'strand', 'start', 'end', 'read_name', ]))

    cram = pysam.AlignmentFile(args.cram, "rc")
    for read in cram.fetch(reference=args.chrom, start=args.start, end=args.end):
        if not read.is_supplementary and not read.is_secondary:
            positions, mod = get_modified_reference_positions(read)
            for pos in positions:
                if pos is not None:
                    print("{chrom}\t{strand}\t{st}\t{en}\t{rname}\t{mod}".format(
                        chrom=read.reference_name,
                        strand='-' if read.is_reverse else '+',
                        st=pos,
                        en=pos,
                        rname=read.query_name,
                        mod=mod))


def get_modified_reference_positions(read):
    basemod = read.get_tag('MM').split(',', 1)[0]
    if '-' in basemod:
        sys.exit("ERROR: modifications on negative strand currently unsupported.")
    base, mod = basemod.split('+')
    deltas = [int(i) for i in read.get_tag('MM').split(',')[1:]]
    locations = np.cumsum(deltas) + np.concatenate((np.zeros(shape=1),
                                                    np.ones(shape=len(deltas) - 1))).astype('int')
    base_index = np.array(
        [i for i, letter in enumerate(read.get_forward_sequence()) if letter == base]
    )
    modified_bases = base_index[locations]
    refpos = np.array(read.get_reference_positions(full_length=True))
    if read.is_reverse:
        return np.flipud(refpos)[modified_bases], basemod
    else:
        return refpos[modified_bases], basemod


def get_args():
    parser = ArgumentParser(description="convert ont-cram file with MM/MP tags to tsv")
    parser.add_argument("cram", help="aligned cram file containing MM/MP tags")
    parser.add_argument("--chrom", help="chromosome to limit calls to.", default=None)
    parser.add_argument("--start", help="start coordinates to limit calls from")
    parser.add_argument("--end", help="end coordinates to limit calls from")
    return parser.parse_args()


if __name__ == '__main__':
    main()
