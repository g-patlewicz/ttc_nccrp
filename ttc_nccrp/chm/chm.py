import numpy as np
from rdkit import Chem
import os
from rdkit.Chem import rdFingerprintGenerator
import pandas as pd
from pathlib import Path
from ttc_nccrp.chm.data import hpc_smarts,dlc_smarts,opc_smarts


class MoleculeProcessor:
    def __init__(self):

        
    def add_chemical(self, dtxsid, smiles):
        
        mol = Chem.MolFromSmiles(smiles)
        self.dtxsid = dtxsid
        if not mol:
            raise ValueError(f'Invalid smiles for {dtxsid}')
        self.molecules[dtxsid]=mol
        return self.molecules
    
    


    def has_metal_atom(self):
        metal_dict = {'Na': 11 , 'Mg': 12, 'Si': 14, 'K':19, 'Ca':20, 'Mn':25, 'Fe':26, 'Cu':29, 'Zn':30, 'Co':27, 'Ni':28, 'As': 33, 'Cr':24, 'Hg':80, 'Pb':82, 'V':23, 'Al':13, 'Ag':47, 'Cd':48, 'B':5, 'Ti': 22, 'Se': 34, 'Sn': 50, 'Sb':51, 'Be': 4, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Te':52, 'Ba':56, 'W':74, 'Au': 79, 'Bi': 83}
        essential_metal_dict = {'Na': 11 , 'K':19,'Mg': 12,  'Ca':20, 'Fe':26,  'Mn':25, 'Co':27, 'Cu':29, 'Zn':30, 'Mo': 42}
        if self.molecules is None:
            return False
        return any(atom for atom in self.molecules.GetAtoms() if atom.GetAtomicNum() in list(metal_dict.values()))

    def metal_ions(self):
        essential_metal_ions = ['[Na+]', '[K+]', '[Mg++]',  '[Ca++]', '[Fe+3]',  '[Mn++]', '[Co+]', '[Cu++]', '[Zn++]', '[Mo++]']
        if self.mol  is None:
            return False
        return any(self.mol.HasSubstructMatch(Chem.MolFromSmarts(e)) for e in essential_metal_ions)

    def P_inorg(self):
        if self.mol  is None:
            return False
        return self.mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]P(=[O])([OH])[OH]'))

    def dlc(self):
        if self.mol is None:
            return False
        return  any([self.mol.HasSubstructMatch(e) for e in dlc_smarts.values()])
    def hpc(self):
        if self.mol is None:
            return False
        return any([self.mol.HasSubstructMatch(e) for e in hpc_smarts.values()])

    def opc(self):
        if self.mol is None:
            return False
        return any([self.mol.HasSubstructMatch(e) for e in opc_smarts.values()])

    def steroid(self):
        if self.mol is None:
            return False
        return self.mol.HasSubstructMatch(Chem.MolFromSmarts('C[C@@]12CCC3c4c(CCC3C1CC[C@H]2O)cc(O)cc4'))


    
        




 