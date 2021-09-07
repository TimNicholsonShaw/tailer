import argparse
import os

try: 
    from Tailer.TailerFunctions import * #hmm, something is weird here, here's a workaround
    from Tailer.PairwiseAligner import *
    # if someone smarter knows what I'm doing wrong, let me know: timnicholsonshaw@gmail.com
except:
    from TailerFunctions import *
    from PairwiseAligner import *

def runGlobal(args):
        # make sql database from gtf, using gene annotations, one time
    db = getOrMakeGTFdb(args.annotation)

    for file in args.files:
        pre, ext = os.path.splitext(file) #Get extension and filename

        print("Tailing file: " + file)

        # convert, if necessary, and index bam
        # no longer necessary to convert to bam and index
        # can refactor this name at some point
        alignments = getHandleOnIndexedBam(file) 

        # convert alignments to a more easily handled dict
        aln_dict = makeAlignmentDict(alignments) 

        print("Calculating tails...")
        # tails everything and creates a dictionary of tailed reads
        tailedReads = makeTailedReadsDict(aln_dict, db)

        # write to file
        tailedReadsToTailFile(tailedReads, pre + "_tail.csv", threeEndThresh=args.threshold)

        print("Wrote " + pre + "_tail.csv " + " to disk.")

def runLocal(args):
    print("Running Locally")
    localAligner(args)

def main():
    """Converts SAM/BAM files to tail files using a gtf"""

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Convert SAM/BAM to a tail file using a gtf")
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-a", "--annotation", help="A GTF formatted annotation (Global Mode)", metavar="")
    group.add_argument("-e", "--ensids", help="Ensembl IDs of genes to query (Local Mode Only)", metavar="")
    parser.add_argument("-t", "--threshold", type=int, default=100, help="Maximum distance from mature end to be included (default=100)", metavar="")
    parser.add_argument("files", nargs="+", help="SAM or BAM formatted files | FASTA/Q for local mode")

    args = parser.parse_args()

    if args.annotation:
        print("Running in global mode")
        runGlobal(args)

    elif args.ensids:
        print("Running in local mode")
        runLocal(args)

    else:
        raise Exception("Something has gone bafflingly wrong with your arguments. This shouldn't be possible. Kudos.")


if __name__ == "__main__":
    main()