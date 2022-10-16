from unicodedata import name
from ase import Atoms, Atom
from ase.io import read, write
from itertools import combinations
from collections import defaultdict
import random
import math
from ase.io.trajectory import Trajectory
from ase.calculators.vasp import Vasp


import os


class MCDFT:

    def __init__(self, atoms, Temp, N, traj) -> None:
        self.atoms=atoms
        self.base_path=os.getcwd()
        self.structures=[atoms]
        self.energy=[]
        self.Temp=Temp
        self.N=N
        self.traj=traj
        self.command='mpirun -n 4 ~/Downloads/vasp/vasp.5.4.4/bin/vasp_std'

    
    def swap_atoms(self, atoms):
        index=defaultdict(list)
        for atom in atoms:
            index[atom.symbol].append(atom.index)
        elements=list(index.keys())
        swap_elements=random.choice(list(combinations(elements, 2)))
        idx1=random.choice(index[swap_elements[0]])
        idx2=random.choice(index[swap_elements[1]])
        atoms[idx1].symbol=swap_elements[1]
        atoms[idx2].symbol=swap_elements[0]
        return atoms

    def monte_carlo(self,dE):
        kB = 1./11604.
        prob = min(math.exp(-dE/(kB*self.Temp)),1.0)
        accept=1
        if prob<random.random():
            accept=0
        return accept
    
    def build_traj(self):
        dir=f"./{self.Temp}k/0"
        atoms=self.structures[-1]
        atoms=self.run_vasp(atoms,dir)
        e=atoms.get_potential_energy()
        self.energy.append(e)

        for i in range(1,self.N):
            dir=f"./{self.Temp}k/{i}"
            atoms=self.swap_atoms(self.structures[-1])
            try:
                atoms=self.run_vasp(atoms,dir)
            except:
                continue
            e=atoms.get_potential_energy()
            dE=e-self.energy[-1]
            accept=self.monte_carlo(dE)
            atoms.info['accept']=accept
            self.traj.write(atoms)
            if accept:
                self.structures.append(atoms)
                self.energy.append(e)

    def run_vasp(self,atoms,dir):
        
        calc = Vasp(directory=dir,
            command=self.command,
                ibrion=-1,
                isif=2,
                ialgo=48,
                nsw=0,
                ismear=0,
                sigma=0.1,
                ediff=0.1e-04,
                prec='Normal',
                xc='PBE',
                encut=420,
                lreal='Auto'
                )
        atoms.pbc=True
        atoms.calc = calc

        atoms.get_potential_energy(atoms)

        return atoms
    
if __name__=="__main__":
    Temp=100
    N=50
    traj=Trajectory(f'TestRun_{Temp}k_{N}steps.traj','w')
    atoms=read('CONTCAR-HEA1-100k')
    mc_run=MCDFT(atoms, Temp, N, traj)
    mc_run.build_traj()
