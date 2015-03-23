import numpy as np
import pandas as pd
import os
import sys
import csv
import pybedtools

"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Splits the 'bedtools coverage' output hist file into two files and computes summary stats for each target

    Parameters
    ----------
    1) Path to the hist file
	2) Path to output summary stats per probe
	3) Path to output summary stats per gene
    Returns
    -------
    
"""


hist_path = sys.argv[1]
probe_stats_path = sys.argv[2]
gene_stats_path = sys.argv[3]
ddg2p_path = sys.argv[4]

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

def gene_coverage(hist_sum_path,gene_bed,output_path):
	histbed = pybedtools.BedTool(hist_sum_path)
	bedbed = pybedtools.BedTool(gene_bed)
	fieldnames = ["symbol","size","cumulative_size_targetted","median_targetted_coverage","nr_targets_below_5x","nr_targets_below_10x","nr_targets_below_20x"]
	csvwriter = csv.DictWriter(open(output_path,'w'), fieldnames=fieldnames,delimiter="\t")
	csvwriter.writeheader()
	for gene in bedbed:
		res = {} 
		hits = histbed.all_hits(gene)
		res["symbol"] = gene.name
		res["size"] = int(gene.end - gene.start)
		res["cumulative_size_targetted"] = 0
		median_bases_read = 0
		res["nr_targets_below_5x"] = 0
		res["nr_targets_below_10x"] = 0
		res["nr_targets_below_20x"] = 0

		for hit in hits:
			hit_length = int(hit.end-hit.start)
			res["cumulative_size_targetted"] += hit_length
			hit_median_coverage = float(hit.name)
			median_bases_read += hit_length*hit_median_coverage
			if hit_median_coverage <= 5:
				res["nr_targets_below_5x"] += 1
			if hit_median_coverage <= 10:
				res["nr_targets_below_10x"] += 1
			if hit_median_coverage <= 20:
				res["nr_targets_below_20x"] += 1

		try:
			res["median_targetted_coverage"] = median_bases_read/res["cumulative_size_targetted"]
		except ZeroDivisionError:
			res["median_targetted_coverage"] = 0 
		csvwriter.writerow(res)

	



if __name__ == "__main__":
	gene_coverage(hist_sum_path,ddg2p_path,output_path)
	(tmp_all_path,tmp_interval_path) = split_hist(sys.argv[1])
	process_intervals(tmp_interval_path,probe_stats_path)
	gene_coverage(probe_stats_path,ddg2p_path,gene_stats_path)
