import unittest
import sys
import os
import tempfile
from pathlib import Path
from Bio import SeqIO
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main import filter_fastq

class TestFileOperations(unittest.TestCase):
    def setUp(self):
        self.test_input = tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fastq',
            delete=False,
            encoding='utf-8'
        )
        test_data = '''@read1
ATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFF
@read2
GCGCGCGCGCGCGCGC
+
HHHHHHHHHHHHHHHH
@read3
AAAAAAAAAAAAAAAA
+
!!!!!!!!!!!!!!!!
@read4
ATATATATATATATAT
+
BBBBBBBBBBBBBBBB'''
        self.test_input.write(test_data)
        self.test_input.close()
        self.test_output = tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fastq',
            delete=False,
            encoding='utf-8'
        )
        self.test_output.close()
    
    def tearDown(self):
        Path(self.test_input.name).unlink(missing_ok=True)
        Path(self.test_output.name).unlink(missing_ok=True)
    
    def test_filtering_by_gc_read_write(self):
        stats = filter_fastq(
            input_file = self.test_input.name,
            output_file = self.test_output.name,
            gc_bounds=(40, 60)
        )
        self.assertEqual(stats['total'], 4)
        self.assertEqual(stats['passed'], 1)
        self.assertTrue(Path(self.test_output.name).exists())
        records = list(SeqIO.parse(self.test_output.name, 'fastq'))
        self.assertEqual(len(records), stats['passed'])
        if len(records) > 0:
            self.assertEqual(str(records[0].seq), 'ATGCATGCATGCATGC')

if __name__ == '__main__':
    unittest.main()