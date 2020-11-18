#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
"""Extract methylation from fast5 files into a sam file. 

Created on Thursday, 25. July 2019.
"""
import numpy as np
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(description="""Extract base modifications from fast5 files called with Guppy 3.3. 


    The modifications are provided as MM and MP tags conforming to https://github.com/samtools/hts-specs/pull/418/commits/11d7fb900b6d51417f59d7cd2cc8540c1b982590
    
    If requested, the modification likelihoods are provided in hex string tags such that ml:Z:C+m,ffff04,A+a,000093. 
        The syntax is 

        ml:Z:([ACGTN][-+][a-z],([0-9a-f][0-9a-f])+);)+

        ml is the custom SAM tag, Z symbol for character string value. Next three groupings are as in 'Base modifications' 
        in  https://github.com/samtools/hts-specs/pull/418 
        Final group is the hexadecimal coded likelihood for the given type modification for each position in SEQ in 
        the original strand. Note that the likelihoods are 'given the underlying base is called the one defined 
        in the tag'. Each hexadecimal value ranging 00-ff is 255*P(data|base, modified).""")
    
    parser.add_argument("input_fast5",nargs="+",
                        help="Input paths of fast5 files  [default:%(default)s]",
                        )
    parser.add_argument("-o", "--output",
                        help="Output the unsorted SAM or passing reads fastq file file here  [default:%(default)s]",
                        default="/dev/stdout")

    parser.add_argument("-f", "--fastq",nargs='?',
                        help="Produce output in fastq format instead of SAM. Store SAM header to file named here. The SAM tags are stored as read comments that can be copied over to SAM [default:%(default)s const:%(const)s]",
                        const="/dev/null")
    parser.add_argument("--failed_reads",
                        help="Output failed reads in fastq format here. With SAM output, the filter/fail is marked with flag  [default:%(default)s]",
                        default="/dev/stderr")                        

    parser.add_argument("-L", "--likelihoods",default=False,action="store_true",
                        help="Include also the raw likelihoods as 'ml' tag. [default:%(default)s]" )

    parser.add_argument("-F", "--filter",default=False,const=7.0, nargs="?",type=float,
                        help="Mark reads with average q less than this as vendor failed [default:%(default)s]" )

    parser.add_argument("-V", "--verbose",default=0,action="count",
                        help="Be more (and more) verbose with output [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose > 0:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')


    logging.info(args)
    return args


class MethylFile(object):


    _KNOWN_MODIFICATIONS = {"5mC":(b'C','m'),
                            "6mA":(b'A','a')}

    def __init__(self,out_strm=None, output_likelihoods=[], filter_mean_q = False, fastq_output=None, failed_reads=None, header=None, verbose=False ):
        """out_strm: Output stream for the sam file, default sys.stdout,
        output_likelihoods:  Output the raw likelihoods for the given modifications"""
        import sys
        import logging as log
        
        self._output_likelihoods = output_likelihoods
        self._filter_mean_q = filter_mean_q
        self.verbose = verbose

        self._fastq_output = fastq_output
        self._fail_strm = failed_reads
        self._header_strm = header


        if out_strm is None:
            self._out_strm = sys.stdout
        else:
            self._out_strm = out_strm
        
        

        self._read_count = 0
        self._quantiles = []
        
        self._cpg_meth = 0  
        self._cpg_total = 0
        
        log.info(self.__dict__)
         
    def write_header(self):
        import json
        import sys
        self._rgid = self._tracking_id["run_id"][-5:]
        self._sample_id = self._tracking_id["sample_id"]
        self._exp_start_time = self._tracking_id["exp_start_time"]
        self._device_id = self._tracking_id["device_id"]
        self._flow_cell_id = self._tracking_id["flow_cell_id"]

        attribs = {"tracking_id":self._tracking_id,
                "basecall_attributes":self._basecall_attributes,
                "modification_attributes":self._mod_attributes}


        self._sam_header = ["@HD\tVN:1.6\tSO:unsorted",
                "@PG\tID:extract_methylation_fast5_to_sam\tPN:extract_methylation_fast5_to_sam.py\tCL:{}""".format(" ".join(sys.argv)),
                "\t".join( [ f"@RG\tID:{self._rgid}",
                        "PL:ONT",
                        f"DT:{self._exp_start_time}",
                        f"PU:{self._flow_cell_id}",f"SM:{self._sample_id}",
                        "on:Z:{}".format(json.dumps(attribs) ) ] )
        ]

        if self._header_strm:
            head_out = self._header_strm
        else:
            head_out = self._out_strm
        head_out.write("\n".join(self._sam_header)+"\n")
        head_out.flush()


    def flush(self):
        self._out_strm.flush()

    def __len__(self):
        "Return number of reads outputted"
        return self._read_count 


    def get_ml_tag(self,mod_base_table):
        """Return modification likelihoods in hex string tags such that ml:Z:C+m,ffff04,A+a,000093. 
        The syntax is 

        ml:Z:([ACGTN][-+][a-z],([0-9a-f][0-9a-f])+);)+

        ml is the custom SAM tag, Z symbol for character string value. Next three groupings are as in 'Base modifications' in  https://github.com/samtools/hts-specs/pull/418
        Final group is the hex coded likelihood for the given type modification for each position in SEQ in the original strand. Note that the likelihoods are
        given when the underlying base is called the one defined in the tag. Each hexadecimal value is 255*P(data|base, modified)
        """
        ml_tag = ["ml:Z:"]
        for base_mod_long in self._output_likelihoods:
            base_mod = self._modified_base_names[base_mod_long]
            mod_likelihoods = mod_base_table[slice(None), base_mod.column_index]
            mod_data = "{}+{},{};".format(base_mod.unmodified_base.decode("ascii"),base_mod.symbol,mod_likelihoods.data.hex())
            ml_tag.append(mod_data)

        return "".join(ml_tag)
 

    def get_mm_mp_tags(self,mod_base_table,SEQ,mod_base_indices={b"C":(3,"m")}):
        "Return list of MM and MP tags conforming to https://github.com/samtools/hts-specs/pull/418/commits/11d7fb900b6d51417f59d7cd2cc8540c1b982590"
        # Require at least phred score 3 for modification to be called. That is at most 66% probability of error.
        # Other option would be phred score 4 with 45% error rate.
        
        import re

        MIN_PHRED = 3 

        tags = []
        MPtag = []
        MMtag = []
        for unmod_base,(mod_index,mod_symbol) in mod_base_indices.items():
            base_indices = np.frombuffer(SEQ.encode("ascii"), dtype=np.uint8)==np.frombuffer(unmod_base,dtype=np.uint8)

            
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


            if self.verbose>1 and unmod_base==b'C' and mod_symbol=="m":
                Mdf = pd.DataFrame(mod_base_table,columns=list(self._mod_attributes["output_alphabet"]))
                Mdf["CALL"] = list(SEQ)
                Mdf["mod_posteriors"] = np.nan
                Mdf.loc[base_indices,"mod_posteriors"] = phred_posteriors
                CpG_idx = [x.start() for x in re.finditer("C(?=G)",SEQ)]
                CpG_posterior = (Mdf.loc[CpG_idx,"mod_posteriors"]).astype(np.uint8)
                #log.info("Likelihoods:CpG:"+str(Mdf.loc[CpG_idx,"Z"]))
                self._cpg_total += len(CpG_idx)
                self._cpg_meth += (CpG_posterior>=MIN_PHRED).sum()

                if (self._read_count+1)%1000 == 0 :
                    log.info("MeanCpGmeth: {}%".format(self._cpg_meth*100.0/self._cpg_total))




                # try to find which likelihood cutoff gives about 70% methylation.
                # Median read will have 70% of CpG sites methylated with Likelihood>=32
                #self._quantiles.append(CpG_likelihood.quantile(0.3))
                #if (self._read_count+1)%1000 == 0 :
                #    log.info("70%meth:"+str(pd.Series(self._quantiles).describe()))

                l_count = CpG_posterior.value_counts()
                if not hasattr(self,"_ll_distrib"):
                    self._ll_distrib = pd.Series(np.zeros(256))

                self._ll_distrib[l_count.index] += l_count
                if (self._read_count+1)%1000 == 0 :
                    X=self._ll_distrib[self._ll_distrib>0].copy()
                    log.info("meth_PhredPosterior_prop:"+(str(X.sort_index().cumsum()/X.sum())))


        if len(MMtag)>0:
            tags.append("MM:Z:"+";".join(MMtag))
            tags.append("MP:Z:"+"".join(MPtag))
        
        return tags

 
    def write_read_fastq(self,QNAME,SEQ,QUAL,tags = [],flags = 0 ):
        if (flags&0x200) == 0:
            outstrm = self._out_strm
        else:
            if self._fail_strm:
                outstrm = self._fail_strm
            else:
                return

        outstrm.write("@"+QNAME)

        for tag in tags:
            outstrm.write("\t")
            outstrm.write(tag)
        outstrm.write("\n")
        outstrm.write(SEQ)
        outstrm.write("\n+\n")
        outstrm.write(QUAL)
        outstrm.write("\n")



    def write_read_sam(self,QNAME,SEQ,QUAL,tags = [],flags = 0 ):
        """mod_base_indices = dict(UNMODIFIED_BASE = (COLUMN_INDEX,BASE_MOD_SYMBOL))
        
        The methylation is reported as described in https://283-3666509-gh.circle-artifacts.com/0/root/project/pdfs/SAMtags.pdf
        """



        # TODO: avoid repeating these
        FLAG=0x4|flags
        RNAME="*"
        POS=0
        MAPQ=255
        CIGAR="*"
        RNEXT="*"
        PNEXT=0
        TLEN=0

        tags = [f"RG:Z:{self._rgid}"] + tags

        line = list(str(x) for x in [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL])

        self._out_strm.write("\t".join(line+tags)+"\n")

    def write_read(self,QNAME,SEQ,QUAL,tags = [],flags = 0 ):
        if self._read_count == 0 :
            self.write_header()

        if self._fastq_output:
            self.write_read_fastq(QNAME,SEQ,QUAL,tags,flags)
        else:
            self.write_read_sam(QNAME,SEQ,QUAL,tags,flags)
        
        self._read_count += 1

    def update_fast5(self,fast5_filepath,analysis_id=None,mod_index=3):
        """Update (i.e. add or change) the methylation data for reads in the given fast5 file.
        
        mod_index gives the index of the modification call table to store in the database. 
                    Default is mC modification. Indices: A,mA,C,mC,G,T = 0,1,2,3,4,5"""
        
        
        from ont_fast5_api.fast5_interface import get_fast5_file
        import numpy as np
        import logging as log
        
        def tqdm(x,*args,**kwargs):
            return x
        if self.verbose>0:
            try:
                from tqdm import tqdm
            except ModuleNotFoundError:
                pass
            

        log.info("Processing file {}".format(fast5_filepath))

        UNMODIFIED_BASES = [b"A", b"A", b"C", b"C", b"G", b"T"]
        assert mod_index >=0 and mod_index < len(UNMODIFIED_BASES), "mod_index must be in the range 0-5."

        BASE = UNMODIFIED_BASES[mod_index]

        log.info("Looking for modification {} of base {}.".format(mod_index,BASE))
        from collections import namedtuple
        BaseMod = namedtuple("BaseMod","column_index shortName longName symbol unmodified_base")
        flags = 0
        with get_fast5_file(fast5_filepath, mode="r") as f5:
            for read_id in tqdm(f5.get_read_ids(),mininterval=30.0):
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
                    self._output_alphabet = read.get_analysis_attributes(analysis_id+"/output_alphabet")

                    _snames = ( (i,x) for i,x in enumerate(self._mod_attributes["output_alphabet"]) if x  not in "ACGT")
                    _lnames =  self._mod_attributes["modified_base_long_names"].split()
                    self._modified_base_names = {lname:BaseMod(i,sname,lname,self._KNOWN_MODIFICATIONS[lname][1], self._KNOWN_MODIFICATIONS[lname][0])  for (i,sname),lname in zip(_snames,_lnames) }

                    log.info(str(self._mod_attributes))
                    log.info(str(self._basecall_attributes))
                    log.info(str(self._tracking_id))
                    log.info(str(self._modified_base_names))


                try:
                    mod_base_table = read.get_analysis_dataset(
                        analysis_id, 'BaseCalled_template/ModBaseProbs')
                except KeyError as err:
                    log.exception(err)
                    log.error("No ModBaseProbs for read {}".format(read_id))
                    continue

                if mod_base_table is None:
                    log.info("No ModBaseProbs for {}".format(read_id))
                    continue

                fastq = read.get_analysis_dataset(
                    analysis_id, 'BaseCalled_template/Fastq')
                if fastq is None:
                    log.info("No Fastq for {}".format(read_id))
                    continue

                try:   
                    seq_title,seq,_,qvals,_ = fastq.split("\n")
                except ValueError as err:
                    log.exception(err)
                    log.error("Bad Fastq for {}:\n{}".format(read_id,fastq))
                    continue

                
                # import pandas as pd
                # Mdf = pd.DataFrame(mod_base_table,index=list(SEQ),columns=list(self._mod_attributes["output_alphabet"]))
                tags = self.get_mm_mp_tags(mod_base_table,seq)

                if len(self._output_likelihoods)>0:
                    tags.append(self.get_ml_tag(mod_base_table))

                if self._filter_mean_q:
                    summary = read.get_summary_data(analysis_id)
                    if summary["basecall_1d_template"]["mean_qscore"] < self._filter_mean_q:
                        flags = 0x200  # QCFAIL   not passing quality controls
                    else:
                        flags = 0 
                
                    

                self.write_read(read_id,seq,qvals,tags,flags)
    
                #assert (self.get(read_id) == mod_likelihoods).all(),"Mismatch on "+read_id

if __name__ == '__main__':
    args=main()

    mdb = MethylFile(open(args.output,"wt"),  output_likelihoods=("5mC","6mA") if args.likelihoods else [],
            filter_mean_q=args.filter,verbose=args.verbose, 
            fastq_output=args.fastq if not args.fastq else open(args.fastq,"wt"),
            header = False if not args.fastq else open(args.fastq,"wt"), 
            failed_reads=None if not args.failed_reads else open(args.failed_reads,"wt"))

    import logging as log
    log.info(args)
    for fn in args.input_fast5:
        mdb.update_fast5(fn)
