import unittest
import os, sys, glob
import gffutils, pysam

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from Tailer import TailerFunctions
    
class TestTailer(unittest.TestCase):
   
    def setUp(self):
        self.alignments = TailerFunctions.getHandleOnIndexedBam(currentdir+"/testSam.sam")
        self.db = TailerFunctions.getOrMakeGTFdb(currentdir+"/test.gtf")
        self.aln_dict = TailerFunctions.makeAlignmentDict(self.alignments)
        self.TailedReads = TailerFunctions.makeTailedReadsDict(self.aln_dict, self.db, rev_comp=False)

    def test_getHandleOnIndexedBam(self):
        #does it create all the files it's supposed to

        self.assertEqual(len(list(self.alignments)), 0) # Generator should be used up by aln_dict

    def test_getOrMakeGTFdb(self):

        # creates DB
        self.assertTrue(os.path.isfile(currentdir +"/test.db"))

        # DB only contains "genes"
        self.assertEqual(len(list(self.db.featuretypes())), 1)
        self.assertEqual(list(self.db.featuretypes())[0], 'gene')

        # DB is the correct length
        #self.assertEqual(len(list(self.db.all_features())), sum(1 for line in open(currentdir + '/test_temp.gtf')))

        # DB is a gffutils object
        self.assertTrue(isinstance(self.db, gffutils.interface.FeatureDB))

        # DB lookup test
        test_gene = self.db["ENSG00000177181"]
        self.assertEqual(test_gene['gene_name'][0], "RIMKLA")
        self.assertEqual(test_gene.start, 42380792)

    def test_attrFinder(self):
        test_line = 'gene_id "ENSG00000160072"; gene_version "20"; gene_name "ATAD3B"; gene_source "ensembl_havana"; gene_biotype "protein_coding";'

        self.assertEqual(TailerFunctions.attrFinder("gene_id", test_line), "ENSG00000160072")
        self.assertEqual(TailerFunctions.attrFinder("gene_name", test_line), "ATAD3B")
        self.assertEqual(TailerFunctions.attrFinder("gene_biotype", test_line), "protein_coding")
    
    def test_makeAlignmentDict(self):

        self.assertEqual(len(self.aln_dict), 4043)
        self.assertTrue(type(self.aln_dict['SRR3643372.856381'][0]), pysam.libcalignedsegment.AlignedSegment)

    def test_getOverlappingGenes(self):

        genes = TailerFunctions.getOverlappingGenes(self.aln_dict["SRR3643372.859002"][0], self.db)
        gene = list(genes)[0]

        self.assertEqual(gene['gene_id'][0], 'ENSG00000230021')
        self.assertTrue(isinstance(gene, gffutils.feature.Feature))

    def test_makeTailedreadsDict(self, rev_comp=False):

        #prove it collapsed reads
        self.assertTrue(len(self.aln_dict) > len(self.TailedReads))
        self.assertTrue(max([self.TailedReads[x].count for x, _ in self.TailedReads.items()])>1)

        # check tail accuracy

    def test_TailClass(self):
        pass

    def test_TailedReadClass(self):
        pass


        # Need some code to test the read_dict object

    def test_tailedReadsToTailFile(self):
        TailerFunctions.tailedReadsToTailFile(self.TailedReads, currentdir+"/out.csv")

        self.assertTrue(os.path.isfile(currentdir+"/out.csv"))

        os.remove(currentdir + "/out.csv")
        
    def tearDown(self):
        for file in glob.glob(currentdir +"/*.bam*"):
            os.remove(file)

        os.remove(currentdir+"/test.db")
class TestLocalEnsID(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass
class TestLocalFASTA(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()