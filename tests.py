import unittest
from noname_tool import count_GGG_triplets


class TestNonameToolMethods(unittest.TestCase):

    def test_count_GGG_triplets(self):
        self.assertEqual(count_GGG_triplets('ATGACCAA'), 0)
        self.assertEqual(count_GGG_triplets(''), 0)
        self.assertEqual(count_GGG_triplets('ATGACCGGGA'), 1)
        self.assertEqual(count_GGG_triplets('AGAGAGGGGAAGGGA'), 3)
        
if __name__ == '__main__':
    unittest.main()
