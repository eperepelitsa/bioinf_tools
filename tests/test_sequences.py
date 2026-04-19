import unittest
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main import DNASequence, RNASequence, AminoAcidSequence

class TestDNASequence(unittest.TestCase):
    def test_dna_complement(self):
        dna = DNASequence('ATGC')
        complement = dna.complement()
        self.assertEqual(str(complement), 'TACG')
        self.assertIsInstance(complement, DNASequence)

    def test_dna_transcription(self):
        dna = DNASequence('ATGC')
        rna = dna.transcribe()
        self.assertEqual(str(rna), 'AUGC')
        self.assertIsInstance(rna, RNASequence)

    def test_dna_reverse_complement(self):
        dna = DNASequence('ATGC')
        rev_comp = dna.reverse_complement()
        self.assertEqual(str(rev_comp), 'GCAT')
        self.assertIsInstance(rev_comp, DNASequence)

class TestRNASequence(unittest.TestCase):
    def test_rna_reverse_transcription(self):
        rna = RNASequence('AUGC')
        dna = rna.reverse_transcribe()
        self.assertEqual(str(dna), 'ATGC')
        self.assertIsInstance(dna, DNASequence)

    def test_rna_complement(self):
        rna = RNASequence('AUGC')
        complement = rna.complement()
        self.assertEqual(str(complement), 'UACG')
        self.assertIsInstance(complement, RNASequence)

class TestAminoAcidSequence(unittest.TestCase):
    def test_peptide_molecular_weight(self):
        peptide = AminoAcidSequence('AAA')
        self.assertAlmostEqual(peptide.calculate_molecular_weight(),
                               231.24, places=2)
        
class TestErrorHandling(unittest.TestCase):
    def test_invalid_dna_seq_error(self):
        with self.assertRaises(ValueError):
            DNASequence('ATGCX')
    
    def test_invalid_rna_seq_error(self):
        with self.assertRaises(ValueError):
            RNASequence('AUGCT')
    
    def test_invalid_protein_seq_error(self):
        with self.assertRaises(ValueError):
            AminoAcidSequence('ABCD')

if __name__ == '__main__':
    unittest.main()