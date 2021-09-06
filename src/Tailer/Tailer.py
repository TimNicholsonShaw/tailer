import argparse
from Tailer.TailerFunctions import * #hmm
import os

def main():
    """Converts SAM/BAM files to tail files using a gtf"""

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Convert SAM/BAM to a tail file using a gtf")

    parser.add_argument("-a", "--annotation", required=True, help="A GTF formatted annotation", metavar="")
    parser.add_argument("-t", "--threshold", type=int, default=100, help="Maximum distance from mature end to be included", metavar="")
    parser.add_argument("files", nargs="+", help="SAM or BAM formatted files")

    args = parser.parse_args()

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


if __name__ == "__main__":
    main()