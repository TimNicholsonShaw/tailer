import os, pysam
from Bio.Seq import Seq
import gffutils
from tqdm import tqdm
import csv

class Tail():
    """Object to hold tail information
    
    Makes use of alignment and gene information to calculate mature end,
    tail length, and tail sequence. Does not save the full alignment or gene
    information for space and speed.
    """
    
    def __init__(self, alignment, gene):
        """Uses alignment and gene objects to calculate tails

        Strips away unimportant information and just leaves the important pieces
        :param alignment    pysam alignment segemnt object
        :param gene         gffutils feature object
        """
        self.readID = alignment.query_name 
        self.is_reverse = alignment.is_reverse # Was the sequence rev comp'd by aligner
        
        # Some genes have no gene name
        try:
            self.gene = gene['gene_name'][0]
        except KeyError:
            self.gene = None
        
        # I don't think these exist, but in case they have no ensembl ID
        try:
            self.geneID = gene["gene_id"][0]
        except KeyError:
            self.geneID = None

        # Calculates mapping in relation to mature end, tail sequence, and tail length
        if alignment.is_reverse: # Positive strand for illumina reads
            read_end = alignment.pos + alignment.reference_length
            self.threeEnd = (read_end) - gene.stop # pysam is zero-indexed, GTFs are 1-indexed

            if alignment.cigartuples[-1][0] == 4: # Interested in softclips at end of cigar, 4=softclipped
                self.tailLen = alignment.cigartuples[-1][1]
                self.tailSeq = alignment.query_sequence[alignment.query_length - self.tailLen:] # saves seq of tail
            else:
                self.tailLen = 0
                self.tailSeq = None
        
        
        else: # Negative strand for illumina reads
            self.threeEnd = gene.start - alignment.pos - 1 # pysam is zero-indexed, GTFs are 1-indexed
            if alignment.cigartuples[0][0] == 4: # Interested in softclips at beginning of CIGAR, 4=softclipped
                self.tailLen = alignment.cigartuples[0][1]
                self.tailSeq = reverse_complement(alignment.query_sequence[:self.tailLen]) #need to rev_comp this because the aligner didn't
            else:
                self.tailLen = 0
                self.tailSeq = None
   
    def __repr__(self):
        return(self.readID + " " 
               + str(self.geneID) + " " 
               + str(self.gene) + " " 
               + str(self.threeEnd) + " " 
               + str(self.tailLen) + " " 
               + str(self.tailSeq))

class TailedRead():
    """Object to hold multiple tails for reads that align multiple places or multiple genes

    Maintains a count of how many times a sequence matching this read has been seen
    And maintains a list of Tail objects for every alignment and gene combination
    """
    
    def __init__(self):
        self.count = 1
        self.tails = []
    
    def __repr__(self):
        return str(self.count) + " reads | " + str(len(self.tails)) + " gene alignments"
    
    def findBestAlignment(self):
        """Returns a list of the "best" tails

        Based on proximity to mature end. Returns the tail object that is closest to the
        zero mature end. In the cases that there are multiple alignments/genes with the
        same mature end, it returns all of them in a list.
        """
        if len(self.tails) == 0: return None

        best = min([abs(x.threeEnd) for x in self.tails]) # Best is tail that's closest to mature end

        out =[]

        for tail in self.tails:
            if abs(tail.threeEnd) == best:
                out.append(tail)
        return out
    

def reverse_complement(seq): 
    """generates a reverse complement of a sequence

    :param seq      string (ATGCN)
    :return         string reverse complemented
    """
    return str(Seq(seq).reverse_complement()) #Used biopython instead

def getHandleOnIndexedBam(samOrBamFile):
    """Reads in a sam or bam file"""
    pre, ext = os.path.splitext(samOrBamFile) #Get extension and filename
    
    # Whether or not it should read the file in binary format
    if ext == ".sam":
        return pysam.AlignmentFile(samOrBamFile, 'r')
    elif ext == ".bam":
        return pysam.AlignmentFile(samOrBamFile, 'rb')

    else:
        raise Exception(ext + " is an unsupported file type")
    
def attrFinder(attrib, attribString):
    """Easy function to pull attributes from column 8 of gtf files"""
    for entry in attribString.split(";"):
        if entry.lstrip().startswith(attrib):
            return entry.lstrip().split(" ")[1][1:-1]
    return None

def getOrMakeGTFdb(GTForDB):
    """Builds or loads a SQL database of gene annotations
    
    :param GTForDB      GTF annotation file or already generated gffutils db
    :return             gffutils database
    """
    #Check if it's a db
    pre, ext = os.path.splitext(GTForDB) #Get extension and filename
    if ext == ".db":
        return gffutils.FeatureDB(GTForDB)
    
    
    #Reduce GTF file to just genes
    with open(GTForDB, "r") as gtffile, open(pre+"_temp.gtf", 'w') as tempfile: #creates a temporary gtf file
        for line in gtffile:
            if line.startswith("#!"): continue #skips comments

            if line.split("\t")[2] == "gene":
                tempfile.write(line)
    
    print("Creating GTF database...")
    db = gffutils.create_db(pre+"_temp.gtf", dbfn=pre+".db", force=True, disable_infer_genes=True, disable_infer_transcripts=True)
    os.remove(pre+"_temp.gtf") # Cleanup temp db, can probably just do this with temporary file library

    return db

def makeAlignmentDict(handledBAM):
    """Organizes alignments by read in a dict"""
    alignment_dict = {}
    
    # places all alignments for a particular read in a list associated with that read name
    for aln in handledBAM:
        alignment_dict[aln.query_name] = alignment_dict.get(aln.query_name, []) + [aln]
    #the handled bam is used up after this
    return alignment_dict

def getOverlappingGenes(aligned_read, gtf_db):
    """Uses gtf database to retrieve all genes that overlap with a read's alignment

    :param aligned_read         pysam aligned segment
    :param gtf_db               previously generated and loaded SQL db
    :return                     list of gffutils feature object genes
    """
    
    # Set strand of read (illumina reads are reversed)
    strand = "-"
    if aligned_read.is_reverse == True: strand="+"
    
    return (gtf_db.region(seqid=aligned_read.reference_name, 
    start=aligned_read.pos, 
    end=aligned_read.pos+aligned_read.reference_length, 
    strand=strand))

def makeTailedReadsDict(alignment_dict, gtf_db):
    """Collapses reads by sequence and creates Tail objects for every alignment/gene

    :param alignment_dict       generated by above function
    :param gtf_db               SQL database generated by above function
    :return                     dictionary with sequence as keys and TailedRead objects as values
    """

    out_dict = {}
    for key, _ in tqdm(alignment_dict.items()):
        
        try: # If the sequence has been seen before, just increments the count
            out_dict[alignment_dict[key][0].query_sequence].count+=1
            continue
        except KeyError: # If it hasn't been seen before, creates a new TailedRead object
            out_dict[alignment_dict[key][0].query_sequence] = TailedRead()
        
        # Iterates over every alignment and gene combination and creates a Tail object
        #This probably leads to some wrong combos, but it gets sorted out later by filters
        for aln in alignment_dict[key]:
            genes = getOverlappingGenes(aln, gtf_db)
            for gene in genes:
                out_dict[alignment_dict[key][0].query_sequence].tails.append(Tail(aln, gene))

    return out_dict

def tailedReadsToTailFile(TailedReads, outLoc, threeEndThresh = 100):
    """Writes TailedReads object to file
    
    :param TailedReads          Previously generated TailedReads object
    :param outLoc               String indictating where it should be written
    :param threeEndThresh       Threshold for inclusion in write, absolute distance to mature end
    """

    out =[]

    for seq, tailedRead in TailedReads.items():
        bestTails = tailedRead.findBestAlignment()

        if not bestTails: continue

        # Joins ensIDs/geneNames together for reads with multiple best tails
        ensIDs = "|".join([str(x.geneID) for x in bestTails])
        geneNames = "|".join([str(x.gene) for x in bestTails])

        # If it's further than this threshold from the mature end, it's probably spurious
        if bestTails:
            if bestTails[0].threeEnd < -threeEndThresh: continue
            elif bestTails[0].threeEnd > threeEndThresh: continue

            if bestTails[0].is_reverse:
                pass
            else:
                seq = reverse_complement(seq) # Need to rev_comp for illumina reads

            out.append([seq, 
            tailedRead.count,
            ensIDs, 
            geneNames, 
            bestTails[0].threeEnd, 
            bestTails[0].tailLen, 
            bestTails[0].tailSeq])
    
    with open(outLoc, 'w') as csvfile:
        writer = csv.writer(csvfile)

        for line in out:
            writer.writerow(line)



    