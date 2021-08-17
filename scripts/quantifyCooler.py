import argparse
import sys
import os

import cooler
import scipy.stats
import numpy as np
import gzip

def opener(filename):
    """
    Determines if a file is compressed or not
    """
    f = open(filename, 'rb')
    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

def to_string(s):
    """
    This convert string or byte to string
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        assert(sys.version_info[0] != 2)
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
    return s

def readPairBEDRegionsToQuantify(bed_pe_file_to_quantify):
    quantif_regions = {}
    with opener(bed_pe_file_to_quantify) as f:
        for line in f:
            line = to_string(line)
            if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                continue
            ls = line.strip().split("\t")
            if len(ls) < 7:
                raise Exception(f"line: {line} does not have 7 fields."
                                " The bedpe file with regions to quantify should be"
                                "chr1 start1 end1 chr2 start2 end2 name.")
            if ls[6] in quantif_regions:
                raise Exception(f"The region {ls[6]} is present twice in {bed_pe_file_to_quantify}.")
            quantif_regions[ls[6]] = (f"{ls[0]}:{ls[1]}-{ls[2]}", f"{ls[3]}:{ls[4]}-{ls[5]}")
    return(quantif_regions)

def quantifyCool(cool1, cool2, quantif_regions, use_raw, fo):
    fo.write("Region\tMean_first_file\tMean_second_file\tratio\tpval_Mann_Whitney_U_test\tpval_Wilcoxon_signed_rank_test\n")
    both_cool_files = [cool1, cool2]
    both_cool = [cooler.Cooler(f) for f in both_cool_files]

    for region in quantif_regions:
        region1, region2 = quantif_regions[region]
        # Data is loaded from cool file
        try:
            data = [c.matrix(balance=(not use_raw)).fetch(region1, region2) for c in both_cool]
        except Exception as e:
            print(f"Could not load the balanced data from {cool1}, {cool2} on {region1}, {region2}")
            raise e

        # When region1 and region2 are identical only the upper matrix is kept
        if region1 == region2:
            data = [d[np.triu_indices(d.shape[0])] for d in data]
        else:
            data = [d.flatten() for d in data]

        # Mean distance is computed
        mean1, mean2 = [np.mean(d) for d in data]
        # U test is performed
        _, u_test = scipy.stats.mannwhitneyu(*data)
        # Wilcoxon
        _, wilcoxon_test = scipy.stats.wilcoxon(*data)
        fo.write(f"{region}\t{mean1}\t{mean2}\t{mean1 / mean2}\t{u_test}\t{wilcoxon_test}\n")

argp = argparse.ArgumentParser(
    description=("Quantify and test statistically pairs of regions on 2 cool files."))
argp.add_argument('--bedpe', default=None, required=True,
                  help="Input a bedpe-like file with pairs of regions to quantify."
                       " Expected format is chr1 start1 end1 chr2 start2 end2 name"
                       " (can be gzip).")
argp.add_argument('--cool1', default=None, required=True,
                  help="First cool file to evaluate.")
argp.add_argument('--cool2', default=None, required=True,
                  help="Second cool file to evaluate.")
argp.add_argument('--raw', action='store_true',
                  help="Use raw values instead of balanced.")
argp.add_argument('--output', default=sys.stdout,
                  type=argparse.FileType('w'),
                  help="Output table with quantification.")

args = argp.parse_args()

quantif_regions = readPairBEDRegionsToQuantify(args.bedpe)
quantifyCool(args.cool1, args.cool2, quantif_regions, args.raw, args.output)
