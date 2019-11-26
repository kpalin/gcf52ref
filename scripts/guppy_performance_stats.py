#!/usr/bin/env python
#-*- coding: utf-8 -*-
 
""" 

Created on Thursday, 21. November 2019.
"""


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate read call performance from guppy")
    
    parser.add_argument("input",nargs="+",
                        help="Input guppy log file [default:%(default)s]",
                        )
    
    parser.add_argument("-V", "--verbose",default=False,action="store_true",
                        help="Be more verbose with output [and log to a file] [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')


    return args

class GuppyLog(object):
    def __init__(self,infile):
        self.infile = infile
        import logging as log
        import re

        self.start_times = dict()
        self.complete_times = dict()

        from datetime import datetime
        p = re.compile("^([^[]+) \[guppy/info\] Read '([0-9a-f-]+)' from file \"([^\"]+)\" (has been loaded|completed).")
        with open(infile,"rt") as f:
            for line in f:
                m = p.match(line)
                if m is None:
                    log.info(line.strip())
                else:
                    t,read_id,f5_name,read_state = m.groups()
                    t = datetime.fromisoformat(t)
                    
                    if read_state ==  'has been loaded':
                        self.start_times[read_id] = t
                    elif read_state =="completed":
                        self.complete_times[read_id] = t
                    else:
                        raise ValueError(m.groups())

        
        self.read_durations = {}
        
        self.common_reads = set(self.start_times.keys())&set(self.complete_times)
        for rname in self.common_reads:
            dur = self.complete_times[rname] -  self.start_times[rname]
            self.read_durations[rname] = dur

    def timespan(self):
        return max(self.complete_times[x] for x in self.common_reads)   - min(self.start_times[x] for x in self.common_reads) 

    def __len__(self):
        return len(self.read_durations)


    def mean_rateHz(self):
        import pandas as pd
        return 1.0/pd.Series(self.read_durations).median().total_seconds()

    def __str__(self):
        import pandas as pd
        return str(pd.Series(self.read_durations).describe())


if __name__ == '__main__':
    args=main()
    for iname in args.input:
        l = GuppyLog(iname)
        if len(l.complete_times)>0:
            print("{}: {:g} reads per second. Total {} reads in {} seconds, i.e. {:g} reads per second average.".format(iname, l.mean_rateHz(),
                len(l),l.timespan().total_seconds(), len(l)/l.timespan().total_seconds()))
            #print(l)