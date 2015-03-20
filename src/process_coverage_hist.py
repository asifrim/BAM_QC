import numpy as np
import pandas as pd
import os
import sys

"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Splits the 'bedtools coverage' output hist file into two files and computes summary stats for each target

    Parameters
    ----------
    1) Path to the hist file
	2) Path to output summary stats per probe
	3) 
    Returns
    -------
    
"""


hist_path = sys.argv[1]
output_path = sys.argv[2]

def split_hist(hist_path):
	tmp_all_path = os.path.dirname(hist_path)+"/"+os.path.basename(hist_path).strip(".gz")+".all.gz"
	tmp_interval_path = os.path.dirname(hist_path)+"/"+os.path.basename(hist_path).strip(".gz")+".intervals.gz"
	command_intervals = "zcat "+hist_path+" | awk  '$1 != \"all\"' | sort | uniq | gzip > "+tmp_interval_path
	command_all = "zcat "+hist_path+" | awk  '$1 == \"all\"' | sort | uniq | gzip > "+tmp_all_path	
	os.popen(command_intervals).read()
	os.popen(command_all).read()
	return (tmp_all_path,tmp_interval_path)

def process_intervals(interval_path,output_path):
	hist = pd.read_table(interval_path,compression="gzip",header=None,names=("chrom","start","stop","depth","bases_at_depth","length","percentage_at_depth","cumsum"),low_memory=False)
	groups = hist.groupby(["chrom","start","stop"])
	output_file = open(output_path,'w')
	for name,group in groups:
		res = {}
		(res["chrom"],res["start"],res["end"]) =  name
		# res["less_than_5x"] = np.sum(group.ix[group.depth <= 5,"percentage_at_depth"])
		# res["higher_than_20x"] = np.sum(group.ix[group.depth >= 20,"percentage_at_depth"])
		group = group.sort("depth")
		group["cumsum"] = np.cumsum(group.ix[:,"percentage_at_depth"])
		res["median"] = group[group["cumsum"] >= 0.5].iloc[0]["depth"]
		
		output_file.write("%(chrom)s\t%(start)s\t%(end)s\t%(median)s\n" % res)

def process_all(all_path):
	hist = pd.read_table(interval_path,compression="gzip",header=None,names=("all","depth","bases_at_depth","length","percentage_at_depth","cumsum"),low_memory=False)
	hist = hist.sort("depth")
	hist["cumsum"] = np.cumsum(hist.ix[:,"percentage_at_depth"])
	res["median"] = hist[hist["cumsum"] >= 0.5].iloc[0]["depth"]	
	print "all\t%(median)s" % res

def gene_coverage(hist_path,gene_interval):
	pass


if __name__ == "__main__":
	(tmp_all_path,tmp_interval_path) = split_hist(sys.argv[1])
	process_intervals(tmp_interval_path,output_path) 