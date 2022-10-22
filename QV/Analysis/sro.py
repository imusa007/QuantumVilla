from ase.io import read, write
from ase.visualize import view
from ase.neighborlist import neighbor_list
from scipy.signal import find_peaks
import numpy as np
from ase.geometry.analysis import Analysis
from collections import Counter
from collections import defaultdict

class SRO:

    def __init__(self, atoms,cutoffs=None):
        
        self.atoms=atoms
        self.cutoffs=cutoffs
        if self.cutoffs is None:
            self.get_shell_cutoff()
        self.cutoffs.sort(reverse=True)
        #print(self.cutoffs)
        self.elements=list(self.atoms.symbols)
        self.atom_dict={k:v for k,v in enumerate(self.elements)}
        self.ndatas=self.get_neighbor_count()
        

    def get_shell_cutoff(self):
        a,b,c,_,_,_=self.atoms.cell.cellpar()
        rmax=min(a,b,c)/2.0-0.5
        ana=Analysis(self.atoms)
        rdfs,distances=ana.get_rdf(rmax,nbins=100,return_dists=True)[0]
        peak_idx=find_peaks(rdfs)[0]
        peak_pos=np.array(distances[peak_idx])
        self.cutoffs=list(((peak_pos[:-1]+peak_pos[1:])/2.0)[:3])

    def get_neighbor_count(self):
        nth_shell_counts=[]
        for cutoff in self.cutoffs:
            i,j=neighbor_list('ij',self.atoms,cutoff)
            data=defaultdict(list)
            for p,q in zip(i,j):
                data[self.atom_dict[p]].append(self.atom_dict[q])
            nth_shell_counts.append({key:Counter(data[key]) for key in data.keys()})
        
        #print(len(nth_shell_counts))

        for idx in range(len(nth_shell_counts)-1):
            for key in nth_shell_counts[idx].keys():
                nth_shell_counts[idx][key].subtract(nth_shell_counts[idx+1][key])
        return nth_shell_counts
    

    def get_sro(self,A,B):
        sros=[]
        for shell_number in range(len(self.cutoffs)): 
            R_A_random=Counter(self.elements)[A]/sum(Counter(self.elements).values())
            #R_B_random=Counter(self.elements)[B]/sum(Counter(self.elements).values())
            R_A_mc=self.ndatas[shell_number][A][B]/sum(self.ndatas[shell_number][A].values())
            sros.append((1-R_A_mc/R_A_random))
        return sros
        