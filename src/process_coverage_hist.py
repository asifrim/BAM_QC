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
	4) Path to output overall stats
	5) Path to DDG2P bed file
    Returns
    -------
    Outputs 3 different files at different levels of granularity: per captured interval, captured intervals grouped by gene, overall
    
"""




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

def process_all(all_path,output_path):
	hist = pd.read_table(all_path,compression="gzip",header=None,names=("all","depth","bases_at_depth","length","percentage_at_depth","cumsum"),low_memory=False)
	hist = hist.sort("depth")
	hist["cumsum"] = np.cumsum(hist.ix[:,"percentage_at_depth"])
	res = {}
	res["median"] = hist[hist["cumsum"] >= 0.5].iloc[0]["depth"]	
	res["0.1"] = hist[hist["cumsum"] >= 0.1].iloc[0]["depth"]
	res["0.9"] = hist[hist["cumsum"] >= 0.9].iloc[0]["depth"]
	of = open(output_path,'w')
	of.write("%s\t%s\t%s\n"%(res["0.1"],res["median"],res["0.9"]))
	of.close

def gene_coverage(hist_sum_path,gene_bed,output_path):
	histbed = pybedtools.BedTool(hist_sum_path)
	bedbed = pybedtools.BedTool(gene_bed)
	ibed = bedbed.intersect(histbed,wo=True)
	print ibed.fn
	d = pd.read_table(ibed.fn, header=None,names=("chrom","start","stop","symbol","chrom2","start2","stop2","median_dp","overlap"),low_memory=False)
	groups = d.groupby(["symbol","chrom","start","stop"])
	
	fieldnames = ["symbol","chr","start","end","size","cumulative_size_targetted","median_targetted_coverage","nr_targets_below_5x","nr_targets_below_10x","nr_targets_below_20x","total_nr_targets"]
	csvwriter = csv.DictWriter(open(output_path,'w'), fieldnames=fieldnames,delimiter="\t")
	csvwriter.writeheader()
	
	for name,group in groups:
		res = {} 
		
		(res["symbol"],res["chr"],res["start"],res["end"]) = name
		res["size"] = int(res["end"] - res["start"])
		res["cumulative_size_targetted"] = 0
		median_bases_read = 0
		res["nr_targets_below_5x"] = 0
		res["nr_targets_below_10x"] = 0
		res["nr_targets_below_20x"] = 0
		res["total_nr_targets"] = 0

		for index, hit in group.iterrows():
			hit_length = hit["stop2"] - hit["start2"]
			res["cumulative_size_targetted"] += hit_length
			median_bases_read += hit_length*hit["median_dp"]
			if hit["median_dp"] <= 5:
				res["nr_targets_below_5x"] += 1
			if hit["median_dp"] <= 10:
				res["nr_targets_below_10x"] += 1
			if hit["median_dp"] <= 20:
				res["nr_targets_below_20x"] += 1
			res["total_nr_targets"] += 1
		try:
			res["median_targetted_coverage"] = median_bases_read/res["cumulative_size_targetted"]
		except ZeroDivisionError:
			res["median_targetted_coverage"] = 0
		csvwriter.writerow(res)

	



if __name__ == "__main__":
	hist_path = sys.argv[1]
	probe_stats_path = sys.argv[2]
	gene_stats_path = sys.argv[3]
	all_stats_path = sys.argv[4]
	ddg2p_path = sys.argv[5]
	(tmp_all_path,tmp_interval_path) = split_hist(sys.argv[1])
	process_all(tmp_all_path,all_stats_path)
	process_intervals(tmp_interval_path,probe_stats_path)
	gene_coverage(probe_stats_path,ddg2p_path,gene_stats_path)
