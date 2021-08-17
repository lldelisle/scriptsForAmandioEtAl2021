import sys
import argparse
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)


class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """
    def __init__(self, expr, msg):
        self.expr = expr
        self.msg = msg

    def __str__(self):
        return self.expr + "\n" + self.msg


def getViewpoints(chrom, fragmentFile, colC, colS, colE, colI,
                  headerSize, viewpointCoo,
                  viewpointRegionInKb, viewpointRegionInFragment):
    viewpoints = []
    if viewpointRegionInFragment is not None:
        previousFrags = [None] * int(viewpointRegionInFragment)
        remaining = -1
    sys.stderr.write("Reading the fragment file.\n")
    with open(fragmentFile, 'r') as f:
        for i, line in enumerate(f):
            if i < headerSize:
                continue
            v = line.split()
            if len(v) < max(colC, colS, colE, colI):
                raise InputError("The ids of columns specified in the input"
                                 " for the fragment file are incompatible with"
                                 f" the line :{line}")
            # Only look at cis:
            if v[colC] != chrom:
                continue
            if viewpointRegionInFragment is None:
                if int(v[colE]) >= viewpointCoo - \
                   int(viewpointRegionInKb) * 1000 \
                   and int(v[colS]) <= \
                   viewpointCoo + int(viewpointRegionInKb) * 1000:
                    viewpoints.append(f"{chrom}_{v[colI]}")
            else:
                if int(v[colE]) >= viewpointCoo and \
                        int(v[colS]) <= viewpointCoo:
                    viewpoints = [fragID for fragID in previousFrags
                                  if fragID is not None]
                    viewpoints.append(f"{chrom}_{v[colI]}")
                    remaining = int(viewpointRegionInFragment)
                elif remaining > 0:
                    viewpoints.append(f"{chrom}_{v[colI]}")
                    remaining -= 1
                elif remaining == -1:
                    previousFrags.pop(0)
                    previousFrags.append(f"{chrom}_{v[colI]}")
    return(viewpoints)


def countCoverageFromViewpoints(validPair, colc1, colc2, colf1, colf2,
                                headerSize, viewpoints):
    """Return a dictionary with the count of occurence of each value
       in col1 or col2 (0-based index)
       when col2 or col1 in viewpoints"""
    outdic = {}
    with open(validPair, 'r') as f:
        for i, line in enumerate(f):
            if i < headerSize:
                continue
            if i % 10000 == 0:
                sys.stderr.write(f"\r{i} lines of valid pairs file analysed.")
            v = line.split()
            if len(v) < max(colc1, colc2, colf1, colf2):
                raise InputError("len(v) < (max(colc1, colc2, colf1, colf2))",
                                 f"{colc1}, {colc2}, {colf1} and {colf2} are "
                                 "incompatible with the number of columns of "
                                 f"the line : {line}")
            val1 = v[colc1] + "_" + v[colf1]
            val2 = v[colc2] + "_" + v[colf2]
            if val1 in viewpoints and val2 not in viewpoints:
                if val2 not in outdic:
                    outdic[val2] = 1
                else:
                    outdic[val2] += 1
            if val2 in viewpoints and val1 not in viewpoints:
                if val1 not in outdic:
                    outdic[val1] = 1
                else:
                    outdic[val1] += 1
    return(outdic)


def fragToBinInCis(chrom, fragmentFile, colC, colS, colE, colI,
                   headerSize, counts, binsize):
    """Convert a coverage in fragments to
       a coverage in fixed binsize
       using proportion when a fragment overlap multiple bins.
       Output an array where
       v[i] is the coverage from binsize * i to binsize * (i + 1)"""
    current_coverage = 0
    coverage = []
    current_bin = 0
    with open(fragmentFile, 'r') as f:
        for i, line in enumerate(f):
            if i < headerSize:
                continue
            if i % 10000 == 0:
                sys.stderr.write(f"\r{i} lines of fragment file analysed.")
            v = line.split()
            if len(v) < max(colC, colS, colE, colI):
                raise InputError("The ids of columns specified in the input"
                                 " for the fragment file are incompatible with"
                                 f" the line :{line}")
            # Only look at cis:
            if v[colC] != chrom:
                continue
            # Fill the previous bins:
            while int(v[colS]) > (current_bin + 1) * binsize:
                coverage.append(current_coverage)
                current_coverage = 0
                current_bin += 1
            # Now fragment start is above i * binsize
            fragID = v[colC] + "_" + v[colI]
            if fragID in counts:
                count = counts[fragID]
                if int(v[colE]) <= (current_bin + 1) * binsize:
                    # The whole fragment is within the same bin
                    current_coverage += count
                else:
                    # It span multiple bins
                    # Current overlap with the i-bin:
                    current_start = int(v[colS])
                    current_end = (current_bin + 1) * binsize
                    # While the end of the bin is included in the fragment:
                    while int(v[colE]) > current_end:
                        # The proportion of the fragment covering the bin
                        # is applied to the count
                        norm_count = count * (current_end - current_start) / \
                            (int(v[colE]) - int(v[colS]))
                        current_coverage += norm_count
                        # And the current coverage is added:
                        coverage.append(current_coverage)
                        current_coverage = 0
                        # The overlap is updated:
                        current_bin += 1
                        current_start = current_bin * binsize
                        current_end = (current_bin + 1) * binsize
                    # When leaving the loop, only part of the bin is covered by
                    # the fragment:
                    current_end = int(v[colE])
                    norm_count = count * (current_end - current_start) / \
                        (int(v[colE]) - int(v[colS]))
                    current_coverage += norm_count
    # Add the last coverage:
    if current_coverage != 0:
        coverage.append(current_coverage)

    return(coverage)


def writeBedGraph(chrom, binsize, coverage, fo):
    """Write in fo a bedgraph with the coordinates corresponding to the binsize
       and the coverage in coverage"""
    for i, value in enumerate(coverage):
        fo.write(f"{chrom}\t{i * binsize}\t{(i + 1) * binsize}\t{value}\n")


def smooth(array, bin):
    lists = [array[i:-(bin - i - 1)]
             for i in range(bin - 1)] + [array[(bin - 1):]]
    middle_mean = [sum(x) / bin for x in zip(*lists)]
    beginning = [sum(array[:(i+1)]) / (i+1) for i in range(bin - 1)]
    end = [sum(array[-(i+1):]) / (i+1) for i in range(bin - 1)]
    end.reverse
    return(beginning + middle_mean + end)


def get_scaling_factor(array, viewpoint_index, exclusion_bin_size):
    around_viewpoint_signal = \
        sum(array[(viewpoint_index - exclusion_bin_size):(viewpoint_index
            + exclusion_bin_size)])
    return((sum(array) - around_viewpoint_signal) / 1000)


def main(args):
    if args.addHeader is not None:
        args.output.write(args.addHeader)
        args.output.write("\n")

    sys.stderr.write("Getting fragment ids for viewpoint.\n")
    try:
        viewpoint_chromosome, viewpoint_pos = \
            args.viewpointCoo.strip().split(':')
    except ValueError:
        raise InputError("viewpoint_chromosome, viewpoint_pos = args"
                         ".viewpointCoo.strip().split(':')",
                         "Wrong format of --viewpointCoo was expecting"
                         " chromosomeName:position")
    try:
        viewpoint_pos = int(viewpoint_pos)
    except ValueError:
        raise InputError("int(viewpoint_pos)",
                         "Wrong format of --viewpointCoo was expecting"
                         " chromosomeName:position."
                         " But conversion of position to integer failed.")

    viewpointList = getViewpoints(viewpoint_chromosome, args.fragmentFile,
                                  args.colForChr - 1, args.colForStart - 1,
                                  args.colForEnd - 1, args.colForID - 1,
                                  args.lineToSkipInFragmentFile,
                                  viewpoint_pos,
                                  args.viewpointRegionInKb,
                                  args.viewpointRegionInFragment)

    sys.stderr.write(f"Got:{viewpointList}\n")
    viewpoint_index = int(viewpoint_pos / (args.binSize * 1000))
    # Compute the coverage per fragment
    sys.stderr.write("Compute coverage from valid pairs.\n")
    coverage_frag = countCoverageFromViewpoints(args.validPair,
                                                args.colChr1 - 1,
                                                args.colChr2 - 1,
                                                args.colFrag1 - 1,
                                                args.colFrag2 - 1,
                                                args.lineToSkipInValidPair,
                                                viewpointList)
    sys.stderr.write("\nBin the coverage.\n")
    # Bin it
    coverage = fragToBinInCis(viewpoint_chromosome, args.fragmentFile,
                              args.colForChr - 1, args.colForStart - 1,
                              args.colForEnd - 1, args.colForID - 1,
                              args.lineToSkipInFragmentFile, coverage_frag,
                              args.binSize * 1000)
    sys.stderr.write("\nGet the scaling factor.\n")
    # Get the scaling_factor
    scaling_factor = get_scaling_factor(coverage, viewpoint_index,
                                        args.exclusionBin)
    sys.stderr.write("Smooth\n")
    # Smooth
    smoothed_coverage = smooth(coverage, args.smoothBins)
    sys.stderr.write("Write output\n")
    # Write the smoothed normalized
    writeBedGraph(viewpoint_chromosome, args.binSize * 1000,
                  [value / scaling_factor for value in smoothed_coverage],
                  args.output)


argp = argparse.ArgumentParser(description="Report a bedgraph of a virtual"
                               " Capture-C from a validPair file and a "
                               "fragment file")
argp.add_argument('--validPair', default=None, required=True,
                  help="A valid pair file containing the fragment ids.")
argp.add_argument('--colChr1', default=None, type=int, required=True,
                  help="The number of the column for the chromosome of read1 "
                  "in the valid pair file.")
argp.add_argument('--colChr2', default=None, type=int, required=True,
                  help="The number of the column for the chromosome of read2 "
                  "in the valid pair file.")
argp.add_argument('--colFrag1', default=None, type=int, required=True,
                  help="The number of the column for the fragment of read1 in "
                  "the valid pair file.")
argp.add_argument('--colFrag2', default=None, type=int, required=True,
                  help="The number of the column for the fragment of read2 in "
                  "the valid pair file.")
argp.add_argument('--lineToSkipInValidPair', default=0, type=int,
                  help="The number line to skip in the valid pair file"
                  " (default 0).")
argp.add_argument('--fragmentFile', default=None, required=True,
                  help="A file containing the coordinates "
                  "of each fragment id.")
argp.add_argument('--colForChr', default=1, type=int,
                  help="The number of the column for the chromosome in the"
                  " fragment file (default 1).")
argp.add_argument('--colForStart', default=2, type=int,
                  help="The number of the column for the start position in the"
                  " fragment file (default 2).")
argp.add_argument('--colForEnd', default=3, type=int,
                  help="The number of the column for the end position in the"
                  " fragment file (default 3).")
argp.add_argument('--colForID', default=4, type=int,
                  help="The number of the column for the fragment id in the"
                  " fragment file (default 4).")
argp.add_argument('--lineToSkipInFragmentFile', default=0, type=int,
                  help="The number line to skip in the fragment file"
                  " (default 0).")
argp.add_argument('--viewpointCoo', default=None, required=True,
                  help="Viewpoint coordinate as chr11:112782000")
argp.add_argument('--binSize', default=1, type=float,
                  help="bin size in kb of the profile (default 1)")
argp.add_argument('--smoothBins', default=5, type=int,
                  help="number of bins used to smooth (default 5)")
argp.add_argument('--exclusionBin', default=5, type=int,
                  help="number of bin to exclude from normalization"
                  " around viewpoint (default 5 means +-5 bins)")
argp.add_argument('--addHeader', default=None,
                  help="Write this line at the begining of the output file.")
argp.add_argument('--output', default=sys.stdout, type=argparse.FileType('w'))

group = argp.add_mutually_exclusive_group()
group.add_argument('--viewpointRegionInKb', default=2.5, type=float,
                   help="Which region to consider as a viewpoint"
                   " around the viewpoint coordinate "
                   "(default 2.5 means +-2.5kb)")
group.add_argument('--viewpointRegionInFragment', default=None, type=int,
                   help="Which region to consider as a viewpoint"
                   " around the viewpoint coordinate (0 means only the "
                   "fragment of the viewpoint, 1 means 1 before, 1 is the"
                   " viewpoint, 1 after)")

args = argp.parse_args()

main(args)
