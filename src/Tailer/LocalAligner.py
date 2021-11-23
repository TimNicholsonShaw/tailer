# https://www.tutorialspoint.com/biopython/biopython_overview_of_blast.htm

from Bio.Blast import NCBIXML
from Bio import SeqIO

from requests.api import get # Used to parse XML output from blast
import subprocess, os, sys, time, requests, json
from tqdm import tqdm
 
import tempfile, csv #I'm doing something wrong here
try:
    import Tailer.TailerFunctions as tf
except:
    import TailerFunctions as tf

class LocalTail:
    def __init__(self, threePrime=None, tailLen=None, tailSeq=None, gene=None, score=None):
        self.threePrime = threePrime
        self.tailLen = tailLen
        self.tailSeq = tailSeq
        self.gene = gene
        self.score = score
    
class LocalTailedRead:

    def __init__(self, seq, count):#, threePrime=None, tailLen=None, tailSeq=None, notes="", gene=""):
        self.seq = seq # read sequence
        self.count = count # how many of them we've seen
        #self.threePrime = threePrime # where the 3' end maps
        #self.tailLen = tailLen # How long the tail is
        #self.tailSeq = tailSeq # what is the sequence of that tail
        #self.notes = notes
        #self.gene = gene

        # Has the read been tailed
        #if threePrime==None: self.tailed = False
        #else: self.tailed = True
        self.potential_tails=[]

    def rev_comp(self):
        self.seq = tf.reverse_complement(self.seq)

    def trim(self, n: int):
        self.seq = self.seq[:-n]

    def addAlignment(self, localTail:LocalTail):
        self.potential_tails.append(localTail)
    def pickBestAlignment(self):
        if len(self.potential_tails) == 0: return None

        best = max([x.score for x in self.potential_tails]) # Best are the ones with the best blast score

        out =[]

        for tail in self.potential_tails:
            if tail.score == best:
                out.append(tail)
        return out   

def parseFASTQ(fastq):
    """
    parses a fastq file into a list of TailedReads
    """

    temp_dict = {}

    with open(fastq, 'r') as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq")):
            seq = str(record.seq)
            temp_dict[seq] = temp_dict.get(seq, 0) + 1

    reads = []
    for key, value in tqdm(temp_dict.items()):
        reads.append(LocalTailedRead(key, count=value))

    return reads

def buildDBFromEID(EID_list, tempfile="temp.fasta"):
    """
    Take a list of EnsIDs and create a blast database
    """
    fasta_dict = getEnsemblSeqs(EID_list)
    if not fasta_dict: raise("Couldn't contact Ensembl")

    with open(tempfile, 'w') as out_file:
        for key, value in fasta_dict.items():
            out_file.write(">" + key + "\n")
            out_file.write(value + "\n")

    subprocess.call(["makeblastdb", "-in", tempfile, "-dbtype", "nucl"])

def queryFormatter(reads, tempfile="temp_query.fasta", rev_comp=False, trim=0):
    """
    create a fasta formatted temporary file
    """
    with open(tempfile, 'w') as handle:
        n=0
        for read in reads:
            if rev_comp: read.rev_comp() # Think I need to add this for gene-specific stuff, flips seq in place
            if trim: read.trim(trim) # trim off barcode
            handle.write(">" + str(n) + ";" + str(read.count) + "\n")
            handle.write(read.seq + "\n")
            n+=1

def alignBlastDB(query, db, outfile):
    """
    Aligns query to properly formatted blast db
    Outputs a temporary XML file
    """
    subprocess.call(["blastn", "-db", db, "-query", query, "-out", outfile, "-outfmt", "5", '-max_target_seqs', '20'])

    return True

def BlastResultsParser(XML_results, reads, expanded_3prime=50, fullName=False):
    """
    Parses XML output from blastn
    Add data to reads object
    """
    records = list(NCBIXML.parse(open(XML_results)))
    
    for record in tqdm(records):
        if len(record.alignments) == 0: continue # skip if there are no alignments
        for i in range(len(record.alignments)):
            query_len = record.query_length
            idx = int(record.query.split(";")[0])
            sbjct_end = record.alignments[i].hsps[0].sbjct_end
            query_end = record.alignments[i].hsps[0].query_end
            target_len = record.alignments[i].length
            score = record.alignments[i].hsps[0].score

            threePrime = sbjct_end - target_len  + expanded_3prime
            tailLen = query_len - query_end
            tailSeq = reads[idx].seq[query_end:]


            if fullName:
                gene_name = record.alignments[i].title
            else:
                gene_name = record.alignments[i].title[record.alignments[0].title.rfind("|")+3:]

            reads[idx].addAlignment(LocalTail(threePrime=threePrime, tailLen=tailLen, tailSeq=tailSeq, gene=gene_name, score=score))


    os.remove(XML_results)
    
    return reads

def tailbuildr(reads, out_loc, seq_out=False):
    """
    creates a .tail file
    """
    reads = sorted(reads, key=lambda x: x.count, reverse=True)

    with open(out_loc, "w") as csvfile:
        writer = csv.writer(csvfile)
        if seq_out:
            writer.writerow(["Sequence", "Count", "EnsID", "Gene_Name", "End_Position", "Tail_Length", "Tail_Sequence" ])
        else:
            writer.writerow(["Count", "EnsID", "Gene_Name", "End_Position", "Tail_Length", "Tail_Sequence" ])


        for read in reads:
            if read.potential_tails:
                bestTail = read.pickBestAlignment()
                out_name = [x.gene for x in bestTail]
                out_name="|".join(out_name)


                if seq_out:
                    writer.writerow([read.seq, read.count, out_name, out_name, bestTail[0].threePrime+bestTail[0].tailLen, bestTail[0].tailLen, bestTail[0].tailSeq])
                else:
                    writer.writerow([read.count, out_name, out_name, bestTail[0].threePrime+bestTail[0].tailLen, bestTail[0].tailLen, bestTail[0].tailSeq])

def getEnsemblSeqs(ID_list, expand_3prime=50):
  server = "https://rest.ensembl.org"
  ext = "/sequence/id"
  headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

  #convert ID list into json
  ID_list = {"ids": ID_list}
  ID_list = json.dumps(ID_list)

  #request sequences of IDs with expanded 3' end equal to expand_3prime
  r = requests.post(server+ext, headers=headers, data=ID_list, params={"expand_3prime":str(expand_3prime)}) 
  
  if not r.ok: #should add more error handling code here
    r.raise_for_status()
    raise("Couldn't contact ensembl")
    sys.exit()

  out = {}
  
  for item in json.loads(r.content): #converts json to dict object
    out[item['id']] = item['seq']

  return out

def localAligner(args):
    tempDir = tempfile.TemporaryDirectory() #Create temporary directory that will be deleted on exit

    # temporary locations for query and database
    queryFile = tempDir.name + "/query.fasta"
    dbFile = tempDir.name + "db.fa"

    args.eids = args.ensids.split(",")

    # Downloads fasta sequences from ensembl and creates BLASTable database
    buildDBFromEID(args.eids, dbFile) 

    for file in args.files:
        pre, ext = os.path.splitext(file) #Get extension and filename
        reads = parseFASTQ(file)

        queryFormatter(reads, queryFile, rev_comp=args.rev_comp, trim=args.trim) # Formats properly for a BLAST search
        print("Aligning...")
        alignBlastDB(queryFile, dbFile, pre+"_temp.xml")

        print("Parsing...")
        reads = BlastResultsParser(pre+"_temp.xml", reads, expanded_3prime=args.mature+50)
        
        tailbuildr(reads, pre+"_tails.csv")

def localFastaAligner(args):
    tempDir = tempfile.TemporaryDirectory() #Create temporary directory that will be deleted on exit

    # temporary locations for query and database
    queryFile = tempDir.name + "/query.fasta"

    # Build blast db from given reference
    subprocess.call(["makeblastdb", "-in", args.fasta, "-dbtype", "nucl"])

    for file in args.files:
        pre, ext = os.path.splitext(file) #Get extension and filename
        reads = parseFASTQ(file)

        queryFormatter(reads, queryFile, rev_comp=args.rev_comp, trim=args.trim) # Formats properly for a BLAST search

        print("Aligning...")
        alignBlastDB(queryFile, args.fasta, pre+"_temp.xml")
        print("Parsing...")
        reads = BlastResultsParser(pre+"_temp.xml", reads, fullName=True, expanded_3prime=args.mature)

        tailbuildr(reads, pre+"_tails.csv")