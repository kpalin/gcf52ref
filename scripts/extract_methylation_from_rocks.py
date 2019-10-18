#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
""" 

Created on Monday, 29. July 2019.
"""


def main():
    import argparse
    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument("-d", "--mod_data",
                        help="Path to the database containing extracted methylation values [default:%(default)s]",
                        default="base_mods.rocksdb")
    
    parser.add_argument("-a", "--align",
                        help="Aligned bam/cram file with the reads [default:%(default)s]",
                        default="aligned.cram")
    
    parser.add_argument("-r", "--reference",
                        help="Genome fasta  [default:%(default)s]",
                        default="reference.fasta")
    parser.add_argument("-o", "--output",
                        help="Output file for the reference achored methylation calls [default:%(default)s]",
                        default="/dev/stdout")    
    parser.add_argument("--no_header",default=False,action="store_true",
                        help="Do not write header to output [default:%(default)s]" )
    parser.add_argument("-V", "--verbose",default=False,action="store_true",
                        help="Be more verbose with output [and log to a file] [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(process)d:%(message)s')
        #if args.verbose!=True:
        #    log_file_handler = logging.FileHandler(args.verbose)
        #    log_file_handler.setFormatter(logging.getLogger().handlers[0].formatter)
        #    logging.getLogger().addHandler(log_file_handler)

    return args

import collections
Methylcall = collections.namedtuple("Methylcall",["chromosome","strand",  "start",   "end",     "read_name",
                                                "log_lik_ratio","log_lik_methylated","log_lik_unmethylated", 
                                                "num_calling_strands", "num_motifs", "sequence"])

class MethylReporter(object):
    def __init__(self,mod_data,aln_file,genome_fasta):
        from extract_methylation_fast5 import MethylDB
        import pysam
        import os.path
        
        assert os.path.exists(mod_data), f"Database {mod_data} must exist!"
        
        self._db_name = mod_data
        self._ref_name = genome_fasta
        self._aln_name = aln_file
        
        self.db = MethylDB(mod_data)
        self.ref = pysam.FastaFile(genome_fasta)
        self.aln = pysam.AlignmentFile(aln_file,"rb")
        
        
    
    def _get_cpg_sites(self,ref_name,ref_start,ref_end,MOTIF="CG"):
        # Reference sites to call
        import re        
        ref_seq = self.ref.fetch(ref_name,ref_start,ref_end)        
        cpg_sites = [ref_start + i.start() for i in re.finditer(MOTIF,ref_seq)]
        return cpg_sites

    def _get_mod_likelihoods(self,read):
        """Return modification likelihoods of the given read as numpy array for each base of the read.
        The parameter `read` is AlignedSegment from pysam and must contain full sequence of the read,
        i.e. no secondary alignments or hard clips."""
        import numpy as np
        import logging as log
        import mappy

        meth = self.db.get(read.query_name)
                
        if meth is None:
            log.info(f"Couldn't find modification data for {read.query_name}")
            return None
        
        
        
        # The read from mapping is given in the reference orientation but the extraction was done in the read orientation.
        if read.is_reverse:
            read_sequence = mappy.revcomp(read.query_sequence)
            MOD_BASE = b"G"
        else:
            read_sequence = read.query_sequence
            MOD_BASE = b"C"
            
        try:
            assert len(meth) == read_sequence.count("C"), "Unexpected number of methylation observations {} instead of {}.\n Read {}".format(len(meth),read_sequence.count("C"),read.tostring())    
        except AssertionError:
            log.info( "Unexpected number of methylation observations {} instead of {}.\n Read {}".format(len(meth),read_sequence.count("C"),read.tostring())    )
            log.info("Query length: {} len(read sequence): {} Inferred read length: {}".format(read.query_length,len(read_sequence),read.infer_read_length()))
            return None
        
        
        
        
        # Tabulate methylation values to an array
        meth_like = np.zeros(read.query_length, np.uint8)

        read_sequence_arr = np.fromstring(read_sequence,dtype="|S1")
        mod_indices =  np.where(read_sequence_arr==b"C")[0]
        meth_like[mod_indices] = meth  # 255*P(Data | base modified, base called) Likelihood of modification for each base of the called sequence (read_sequence).

        return meth_like
    
    
    
    def _map_calls_to_reference(self,read,string_output=False):
        "Slow implementation of getting the CpG methylation calls from the methylation database"
        import logging as log

        import numpy as np
        import mappy
        
        if read.query_length < read.infer_read_length():
            #log.info(f"Don't have the complete sequence at hand for {read.query_name}")
            return None

        meth_like = self._get_mod_likelihoods(read)
        if meth_like is None:
            log.warning(f"Couldn't find methylation for read {read.query_name}")
            return None
        
        
        cpg_sites = self._get_cpg_sites(read.reference_name,read.reference_start,read.reference_end)

        
        # Map the reference CpG:s to query to read coordinates.
        aln_pairs = read.get_aligned_pairs()

        
        r2q_map = {r_pos:q_pos for q_pos,r_pos in aln_pairs if q_pos is not None }
        
        #log.info("Proportion CpGs aligned {}: {}".format(read.is_reverse,np.array([r in r2q_map for r in cpg_sites]).mean() ))
        

        cpg_sites_in_query = [ (r2q_map[r],r) for r in cpg_sites if r in r2q_map]
        
        
        if read.is_reverse:
            cpg_sites_in_read = [(q+1,r) for (q,r) in cpg_sites_in_query]
            meth_like = meth_like[::-1]
        else:
            cpg_sites_in_read = cpg_sites_in_query
        
 



        
        if string_output:
            fmt = "{}\t{}\t{{}}\t{{}}\t{}\t{{:.3}}\t{{:.3}}\t{{:.3}}\t{}\t{}\t{}\n".format(read.reference_name, "-" if read.is_reverse else "+", 
                                                                                             read.query_name, 1,1,"CG" )
            cpg_sites_in_read_meth = ( (query_pos,ref_pos,np.log10( (meth_like[query_pos]+1.0)/256),np.log10( (256-meth_like[query_pos])/256.0)) for query_pos,ref_pos in  cpg_sites_in_read)
            MOD_EPS = 1.0/256
            ret =  ( fmt.format(ref_pos,ref_pos+1,  
                        log_mod_ll-log_nonmod_ll, 
                        log_mod_ll, 
                        log_nonmod_ll) for query_pos,ref_pos,log_mod_ll,log_nonmod_ll in  cpg_sites_in_read_meth)
            #log.info("Mean methylation for {} {}: {}".format(read.query_name,read.is_reverse, np.array([x.log_lik_methylated for x in ret   ]).mean()))
        else:
            
            methylcall_dtype = [("chromosome",'U5'), ("strand",'U1'),  ("start",np.uint64),   ("end",np.uint64), ("read_name",'U40'),
                                ("log_lik_ratio",np.float32), ("log_lik_methylated",np.float32), ("log_lik_unmethylated",np.float32), 
                                ("num_calling_strands",np.uint8), ("num_motifs",np.uint8), ("sequence",'U3') ]

            ret = np.fromiter( ( (read.reference_name, "-" if read.is_reverse else "+", ref_pos,ref_pos+1, read.query_name,  
                        np.log10((meth_like[query_pos]+1)/(256-meth_like[query_pos])), 
                        np.log10((meth_like[query_pos]+1)/256.0), 
                        np.log10((256-meth_like[query_pos])/256.0),1,1,"CG",meth_like[query_pos]) for query_pos,ref_pos in  cpg_sites_in_read),
                        dtype = methylcall_dtype, count=len(cpg_sites_in_read) ) 
                    #ret = np.fromiter( ( (read.reference_name, "-" if read.is_reverse else "+", ref_pos,ref_pos+1, read.query_name,  
        #            np.log10((meth_like[query_pos]+1)/(256-meth_like[query_pos])), 
        #            np.log10((meth_like[query_pos]+1)/256.0), 
        #            np.log10((256-meth_like[query_pos])/256.0),1,1,"CG",meth_like[query_pos]) for query_pos,ref_pos in  cpg_sites_in_read),
        #            dtype = methylcall_dtype, count=len(cpg_sites_in_read) )       
        #ret = [ Methylcall(read.reference_name, "-" if read.is_reverse else "+", ref_pos,ref_pos+1, read.query_name,  
        #            np.log10((meth_like[query_pos]+1)/(256-meth_like[query_pos])), 
        #            np.log10((meth_like[query_pos]+1)/256.0), 
        #            np.log10((256-meth_like[query_pos])/256.0),1,1,"CG",meth_like[query_pos]) for query_pos,ref_pos in  cpg_sites_in_read]
        
    
        return ret
    
    def __len__(self):
        return len(self.db)
    
    @property
    def manager(self):
        if not hasattr(self,"_manager"):
            import multiprocessing as mp
            self._manager = mp.Manager()
        return self._manager
    
    @staticmethod
    def _read_queue_reads(aln_file,contig,q):
        import pysam
        import logging as log
        log.info(f"Running {aln_file} {contig} {q}")
        
        aln = pysam.AlignmentFile(aln_file,"rb")
        import dill
        for read in aln.fetch(contig):
            if not read.is_secondary:
                q.put(dill.dumps(read))
                
        return contig
            
    
    def _iter_reads(self):
        import multiprocessing as mp
        import logging as log
        from queue import Empty

        q = self.manager.Queue(1000)
        
        refs = list(self.aln.references)
        log.info("Reading reads from reference sequences: "+", ".join(refs))
        import dill
        with mp.Pool(3) as pool:
            _res = pool.starmap_async(self._read_queue_reads,[(self._aln_name,ref_name,q) for ref_name in refs])
            
            while not _res.ready():
                try:
                    yield dill.loads(q.get(timeout=3))
                except Empty:
                    log.info("Stalling.. alignment reading is taking time.")
            log.info(_res.get())
        
    
    
    def fetch_methylation_calls(self,string_output=False):
        import logging as log
        aln_it = iter(self.aln)#.fetch()
        #aln_it = self._iter_reads()
        
        for read in aln_it:
            if read.is_secondary:
                #log.info(f"Skipping secondary alignment of {read.query_name} at {read.reference_name}:{read.reference_start}")
                continue
            yield self._map_calls_to_reference(read,string_output)
            
          
    

if __name__ == '__main__':
    args=main()
    import sys
    import logging as log
    from tqdm import tqdm

    meth_rep = MethylReporter(args.mod_data,args.align,args.reference)
    
    
    
    log.info("Using a database of about {:,} reads.".format(len(meth_rep)))
    n_out = 0
    with open(args.output,"w") as outstrm:
        if not args.no_header:
            outstrm.write("#"+"\t".join(Methylcall._fields)+"\n")
        #fmt = "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\n" 
        fmt = "%s\t%s\t%d\t%d\t%s\t%.3g\t%.3g\t%.3g\t%d\t%d\t%s\t%d"
        import numpy 
        for c in tqdm(meth_rep.fetch_methylation_calls(string_output=True),mininterval=5,unit="reads"):
            if c is None:
                continue            
            #numpy.savetxt(outstrm,c,fmt=fmt)
            #outstrm.writelines( fmt.format(*x) for x in c)
            outstrm.writelines(c)
            n_out +=1
            #if n_out>100000:
            #    break
