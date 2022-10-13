from ase import Atoms, Atom
from ase.io import read, write
from ase.visualize import view
from ase.neighborlist import neighbor_list
from collections import Counter
from collections import defaultdict, OrderedDict
from itertools import combinations, combinations_with_replacement
import random
from ase.build import bulk
import math
from ase.io.trajectory import Trajectory

class SRO:

    def __init__(self, atoms,cutoffs) -> None:
        self.atoms=atoms
        self.cutoffs=cutoffs
        self.elements=list(self.atoms.symbols)
        self.atom_dict={k:v for k,v in enumerate(self.elements)}
        self.ndata=self.get_neighbor_count()

    def get_neighbor_count(self):
        
        i,j=neighbor_list('ij',self.atoms,self.cutoffs)
        data=defaultdict(list)
        for p,q in zip(i,j):
            data[self.atom_dict[p]].append(self.atom_dict[q])
        ndata={key:Counter(data[key]) for key in data.keys()}
        return ndata

    def get_sro(self,A,B): 
        R_A_random=Counter(self.elements)[A]/sum(Counter(self.elements).values())
        R_B_random=Counter(self.elements)[B]/sum(Counter(self.elements).values())
        R_A_mc=self.ndata[A][B]/sum(self.ndata[A].values())
        sro=(1-R_A_mc/R_A_random)
        return sro
        