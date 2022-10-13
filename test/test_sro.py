import unittest
from ase.io import read
from context import QV
from QV.Analysis.sro import SRO
import numpy as np

class TestSRO(unittest.TestCase):

    def test_sro(self):
        atoms=read('CONTCAR-HEA1-T100')
        sro_analyzer=SRO(atoms,cutoffs=2.90)
        
        sro_WW=sro_analyzer.get_sro('W','W')
        sro_TiTa=sro_analyzer.get_sro('Ti','Ta')
        sro_TaCr=sro_analyzer.get_sro('Ta','Cr')
        
        assert np.allclose([sro_WW,sro_TiTa,sro_TaCr], [0.03125, 0.25758, -0.93798], rtol=1e-03, atol=1e-03)

if __name__ == '__main__':
    unittest.main()