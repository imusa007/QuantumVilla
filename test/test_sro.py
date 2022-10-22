import unittest
from ase.io import read
from context import QV
from QV.Analysis.sro import SRO
import numpy as np

class TestSRO(unittest.TestCase):

    def test_sro(self):
        atoms=read('CONTCAR-HEA1-T100')
        sro_analyzer=SRO(atoms) #,cutoffs=[3.60,2.90])
        
        sros_WW=sro_analyzer.get_sro('W','W')
        sros_TiTa=sro_analyzer.get_sro('Ti','Ta')
        sros_TaCr=sro_analyzer.get_sro('Ta','Cr')
        # print(sros_WW[0],sros_TiTa[0],sros_TaCr[0])
        # print(sros_WW[1],sros_TiTa[1],sros_TaCr[1])
        # print(sros_WW[2],sros_TiTa[2],sros_TaCr[2])
        assert np.allclose([sros_WW[-1],sros_TiTa[-1],sros_TaCr[-1]], [0.03125, 0.25758, -0.93798], rtol=0.03, atol=0.03)

if __name__ == '__main__':
    unittest.main()