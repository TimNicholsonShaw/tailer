import argparse
import os
import tempfile
try: 
    import Tailer.TailerFunctions as tf #hmm, something is weird here, here's a workaround
    import Tailer.LocalAligner as la
    # if someone smarter knows what I'm doing wrong, let me know: timnicholsonshaw@gmail.com
except:
    import TailerFunctions as tf
    import LocalAligner as la

def runGlobal(args):
        # make sql database from gtf, using gene annotations, one time
    tempDir = tempfile.TemporaryDirectory() #Create temporary directory that will be deleted on exit

    db = tf.getOrMakeGTFdb(args.annotation)

    for file in args.files:
        pre, ext = os.path.splitext(file) #Get extension and filename

        print("Tailing file: " + file)

        # convert, if necessary, and index bam
        # no longer necessary to convert to bam and index
        # can refactor this name at some point
        alignments = tf.getHandleOnIndexedBam(file) 

        # convert alignments to a more easily handled dict
        aln_dict = tf.makeAlignmentDict(alignments, args.read) 

        print("Calculating tails...")
        # tails everything and creates a dictionary of tailed reads
        tailedReads = tf.makeTailedReadsDict(aln_dict, db, args.rev_comp)

        # write to file
        tf.tailedReadsToTailFile(tailedReads, pre + "_tail.csv", threeEndThresh=args.threshold, seq_out=args.sequence, mature_end=args.mature)

        print("Wrote " + pre + "_tail.csv " + " to disk.")

def main():
    """Handles argument parsing and run selection"""

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Convert SAM/BAM to a tail file using a gtf")
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument("-a", "--annotation", help="A GTF formatted annotation (Global Mode Only), excludes -e", metavar="")
    group.add_argument("-e", "--ensids", help="Ensembl IDs of genes to query, comma separated, no space (Local Mode Only), excludes -a", metavar="")
    group.add_argument("-f", "--fasta", help="FASTA file to act as a reference [Local Only]", metavar="")
    parser.add_argument("-t", "--threshold", type=int, default=100, help="Maximum distance from mature end to be included (default=100)", metavar="")
    parser.add_argument("-r", "--rev_comp", action="store_true",help="Reverse complement reads")
    parser.add_argument("-x", "--trim", default=0, type=int, help="trim off x nucleotides from end [Helper for local only]", metavar="")
    parser.add_argument("-s", "--sequence", action="store_true", help="for debugging, outputs nucleotide sequences to file")
    parser.add_argument("files", nargs="+", help="SAM or BAM formatted files | FASTA/Q for local mode")
    parser.add_argument("-read", "--read", default=2, type=int, help="Which read to use")
    parser.add_argument("-m", "--mature", default=0, type=int, help="Mature end adjustment, most useful for custom fastas in local mode.")


    args = parser.parse_args()

    
    if args.annotation:
        print("Running in global mode")
        runGlobal(args)

    elif args.ensids:
        print("Running in local mode")
        la.localAligner(args)

    elif args.fasta:
        print("Running in local mode with reference FASTA")
        la.localFastaAligner(args)

    else:
        raise Exception("Something has gone bafflingly wrong with your arguments. This shouldn't be possible. Kudos.")


if __name__ == "__main__":
    main()