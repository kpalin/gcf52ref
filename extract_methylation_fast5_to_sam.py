#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
"""Extract methylation from fast5 files into a RocksDB file. Also has an interface to read those values.

Created on Thursday, 25. July 2019.
"""


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract methylation from fast5 files")
    
    parser.add_argument("input_fast5",nargs="+",
                        help="Input paths of fast5 files  [default:%(default)s]",
                        )
    parser.add_argument("-o", "--output",
                        help="Output the unsorted SAM file here  [default:Standard output]",
                        default=None)
    #parser.add_argument("-p", "--processes",type=int,
    #                    help="Database to store the modifications to  [default:%(default)s]",
    #                    default=1)
    parser.add_argument("-V", "--verbose",default=False,action="store_true",
                        help="Be more verbose with output [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')


    return args


class MethylFile(object):
    def __init__(self,out_strm=None ):
        import sys
        if out_strm is None:
            self._out_strm = sys.stdout
        else:
            self._out_strm = out_strm

        self._read_count = 0
        self._quantiles = []
        import pandas as pd
        import numpy as np
        self._ll_distrib = pd.Series(np.zeros(256))
        self._cpg_meth = 0  
        self._cpg_total = 0
         
    def write_header(self):
        import json
        import sys
        self._rgid = self._tracking_id["run_id"][-5:]
        self._sample_id = self._tracking_id["sample_id"]
        self._exp_start_time = self._tracking_id["exp_start_time"]
        self._device_id = self._tracking_id["device_id"]
        self._flow_cell_id = self._tracking_id["flow_cell_id"]
        self._distribution_version = self._tracking_id["distribution_version"] 

        self._sam_header = ["@HD\tVN:1.6\tSO:unsorted",
                "@PG\tID:extract_methylation_fast5_to_sam\tPN:extract_methylation_fast5_to_sam.py\tCL:{}""".format(" ".join(sys.argv)),
                f"@RG\tID:{self._rgid}\tPL:ONT\tDT:{self._exp_start_time}\tPU:{self._flow_cell_id}\tSM:{self._sample_id}",
                "@CO\ttracking_id={}".format(json.dumps(self._tracking_id)),
                "@CO\tbasecall_attrib={}".format(json.dumps(self._basecall_attributes)),
                "@CO\tmodification_attrib={}".format(json.dumps(self._mod_attributes))
        ]

        
        self._out_strm.write("\n".join(self._sam_header)+"\n")


    def flush(self):
        self._out_strm.flush()

    def __len__(self):
        "Return number of reads outputted"
        return self._read_count
    
    def write_read(self,QNAME,SEQ,QUAL,mod_base_table,mod_base_indices={b"C":(3,"m")}):
        """mod_base_indices = dict(UNMODIFIED_BASE = (COLUMN_INDEX,BASE_MOD_SYMBOL))
        
        The methylation is reported as described in https://283-3666509-gh.circle-artifacts.com/0/root/project/pdfs/SAMtags.pdf
        """
        import numpy as np

        if self._read_count == 0 :
            self.write_header()

        # TODO: avoid repeating these
        FLAG=0x4
        RNAME="*"
        POS=0
        MAPQ=255
        CIGAR="*"
        RNEXT="*"
        PNEXT=0
        TLEN=0

        tags = [f"RG:Z:{self._rgid}"]

        line = list(str(x) for x in [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL])
        #end TODO
        
        # Require at least phred score 3 for modification to be called. That is at most 66% probability of error.
        # Other option would be phred score 4 with 45% error rate.

        MIN_PHRED = 3 

        MPtag = []
        MMtag = []
        for unmod_base,(mod_index,mod_symbol) in mod_base_indices.items():
            base_indices = np.fromstring(SEQ,"|S1")==unmod_base

            
            mod_likelihoods = mod_base_table[base_indices, mod_index]
            # mod_likelihoods[i] = 255*P(i:th unmod_base character of SEQ is modified|data)


            # Prior of 70% methylation on CpG sites.
            p_modified = 0.9
            q_unmodified = 1-p_modified
            mod_posteriors = mod_likelihoods*p_modified/ (p_modified*mod_likelihoods+q_unmodified*(256-mod_likelihoods))
            phred_posteriors = np.round(-10*np.log10(1-mod_posteriors)).astype(np.int8)
            #phred_likelihoods = np.round(-10*np.log10((256.0-mod_likelihoods)/256)).astype(np.int8)
            # phred_likelihoods = phred scaled probability of base not being modified.
            
            mod_base_loci = (phred_posteriors >= MIN_PHRED).nonzero()[0]
            if len(mod_base_loci)== 0:
                #log.info(f"No modified bases for {QNAME}")
                continue
            skips = unmod_base.decode("ascii")+f"+{mod_symbol},{mod_base_loci[0]}"

            mod_base_skips = np.diff(mod_base_loci)-1
            if len(mod_base_skips)>0:
                skips += "," +",".join(str(x) for x in mod_base_skips)
            MMtag.append(skips)

            MPtag.extend(chr(x) for x in (phred_posteriors[mod_base_loci]+33))


            if True and unmod_base==b'C' and mod_symbol=="m":
                import pandas as pd
                import re
                Mdf = pd.DataFrame(mod_base_table,columns=list(self._mod_attributes["output_alphabet"]))
                Mdf["CALL"] = list(SEQ)
                Mdf["mod_posteriors"] = np.nan
                Mdf.loc[base_indices,"mod_posteriors"] = phred_posteriors
                CpG_idx = [x.start() for x in re.finditer("C(?=G)",SEQ)]
                CpG_likelihood = (Mdf.loc[CpG_idx,"mod_posteriors"]).astype(np.uint8)
                #log.info("Likelihoods:CpG:"+str(Mdf.loc[CpG_idx,"Z"]))
                self._cpg_total += len(CpG_idx)
                self._cpg_meth += (CpG_likelihood>=MIN_PHRED).sum()

                if (self._read_count+1)%1000 == 0 :
                    log.info("MeanCpGmeth: {}%".format(self._cpg_meth*100.0/self._cpg_total))




                # try to find which likelihood cutoff gives about 70% methylation.
                # Median read will have 70% of CpG sites methylated with Likelihood>=32
                #self._quantiles.append(CpG_likelihood.quantile(0.3))
                #if (self._read_count+1)%1000 == 0 :
                #    log.info("70%meth:"+str(pd.Series(self._quantiles).describe()))

                l_count = CpG_likelihood.value_counts()
                self._ll_distrib[l_count.index] += l_count
                if (self._read_count+1)%1000 == 0 :
                    X=self._ll_distrib[self._ll_distrib>0].copy()
                    log.info("meth_LL_prop:"+(str(X.sort_index().cumsum()/X.sum())))


        if len(MMtag)>0:
            tags.append("MM:Z:"+";".join(MMtag))
            tags.append("MP:Z:"+"".join(MPtag))

        self._out_strm.write("\t".join(line+tags)+"\n")

        
        self._read_count += 1

    def update_fast5(self,fast5_filepath,analysis_id=None,mod_index=3,verbose=False):
        """Update (i.e. add or change) the methylation data for reads in the given fast5 file.
        
        mod_index gives the index of the modification call table to store in the database. 
                    Default is mC modification. Indices: A,mA,C,mC,G,T = 0,1,2,3,4,5"""
        
        
        from ont_fast5_api.fast5_interface import get_fast5_file
        import numpy as np
        import logging as log
        if verbose:
            from tqdm import tqdm
        else:
            def tqdm(x):
                return x

        log.info("Processing file {}".format(fast5_filepath))

        UNMODIFIED_BASES = [b"A", b"A", b"C", b"C", b"G", b"T"]
        assert mod_index >=0 and mod_index < len(UNMODIFIED_BASES), "mod_index must be in the range 0-5."

        BASE = UNMODIFIED_BASES[mod_index]

        log.info("Looking for modification {} of base {}.".format(mod_index,BASE))

        with get_fast5_file(fast5_filepath, mode="r") as f5:
            for read_id in tqdm(f5.get_read_ids()):
                #if read_idx%100:
                #    log.info("Processing read {}".format(read_id))

                read = f5.get_read(read_id)
                if analysis_id == None:
                    analysis_id = read.get_latest_analysis('Basecall_1D')
                    log.info("Using basecall anaysis: {}".format(analysis_id))

                if self._read_count == 0:
                    self._mod_attributes = read.get_analysis_attributes(analysis_id + "/BaseCalled_template/ModBaseProbs")
                    self._basecall_attributes = read.get_analysis_attributes(analysis_id)
                    self._tracking_id = read.get_tracking_id()
                    log.info(str(self._mod_attributes))
                    log.info(str(self._basecall_attributes))
                    log.info(str(self._tracking_id))


                mod_base_table = read.get_analysis_dataset(
                    analysis_id, 'BaseCalled_template/ModBaseProbs')
                if mod_base_table is None:
                    log.info("No ModBaseProbs for {}".format(read_id))
                    continue

                fastq = read.get_analysis_dataset(
                    analysis_id, 'BaseCalled_template/Fastq')
                if fastq is None:
                    log.info("No Fastq for {}".format(read_id))
                    continue
                    
                seq_title,seq,_,qvals,_ = fastq.split("\n") 

                
                # import pandas as pd
                # Mdf = pd.DataFrame(mod_base_table,index=list(SEQ),columns=list(self._mod_attributes["output_alphabet"]))
                self.write_read(read_id,seq,qvals,mod_base_table)
    
                #assert (self.get(read_id) == mod_likelihoods).all(),"Mismatch on "+read_id

if __name__ == '__main__':
    args=main()
    mdb = MethylFile(args.output)

    import logging as log
    log.info(args)
    for fn in args.input_fast5:
        mdb.update_fast5(fn,verbose=args.verbose)
