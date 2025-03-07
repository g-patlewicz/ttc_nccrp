import numpy as np
from rdkit import Chem
import os
import pandas as pd
from pathlib import Path
from ttc_nccrp.chm.data import hpc_smarts,dlc_smarts,opc_smarts, genetox_smarts


class MoleculeProcessor:
    def __init__(self):

        self.molecules = {}
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
        if self.molecules[self.dtxsid] is None:
            return False
        elif any(atom for atom in self.molecules[self.dtxsid].GetAtoms() if atom.GetAtomicNum() in list(metal_dict.values())):
            return 'inorganic'
       

    def metal_ions(self):
        essential_metal_ions = ['[Na+]', '[K+]', '[Mg++]',  '[Ca++]', '[Fe+3]',  '[Mn++]', '[Co+]', '[Cu++]', '[Zn++]', '[Mo++]']
        if self.molecules[self.dtxsid]  is None:
            return False
        elif any(self.molecules[self.dtxsid].HasSubstructMatch(Chem.MolFromSmarts(e)) for e in essential_metal_ions):
            return 'essential metal ion'
        

    def P_inorg(self):
        if self.molecules[self.dtxsid]  is None:
            return False
        elif self.molecules[self.dtxsid].HasSubstructMatch(Chem.MolFromSmarts('[OH]P(=[O])([OH])[OH]')):
            return 'inorganic'
       

    def dlc(self):
        if self.molecules[self.dtxsid] is None:
            return False
        elif any(self.molecules[self.dtxsid].HasSubstructMatch(e) for v in dlc_smarts.values() for e in v):
            return 'Dioxin-like'
        
    def hpc(self):
        if self.molecules[self.dtxsid] is None:
            return False
        elif any(self.molecules[self.dtxsid].HasSubstructMatch(e) for v in hpc_smarts.values() for e in v):
            return 'HPC'
        

    def opc(self):
        if self.molecules[self.dtxsid] is None:
            return False
        elif any(self.molecules[self.dtxsid].HasSubstructMatch(e) for v in opc_smarts.values() for e in v):
            return 'OP or carbamate'
       

    def steroid(self):
        if self.molecules[self.dtxsid] is None:
            return False
        elif self.molecules[self.dtxsid].HasSubstructMatch(Chem.MolFromSmarts('C[C@@]12CCC3c4c(CCC3C1CC[C@H]2O)cc(O)cc4')):
            return 'Steroid'
        

    def genetox_alerts(self):
    # Define the primary pattern for genetox alert
        pattern = 'a[N]=[N]a'
        pattern_mol = Chem.MolFromSmarts(pattern)

        # Define sub-patterns that must NOT be present
        sub_patterns = [
            '[$(a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]N=N[$(a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]',
            '[$(a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]N=N[$(a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]',
            '[$(a:a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]N=N[$(a:a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]',
            '[$(a:a:a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]N=N[$(a:a:a:a:a(S(=[OX1])(=[OX1])([O-,OX2H1])))]'
        ]

        # Retrieve the molecule object safely
        mol = self.molecules.get(self.dtxsid)
        if not mol:
            return False  # No molecule available, return False

        # Check for primary genetox pattern
        if mol.HasSubstructMatch(pattern_mol):
            # Ensure none of the exclusion sub-patterns are present
            if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(sp)) for sp in sub_patterns):
                return 'Genetox'

        # Check additional genetox SMARTS patterns stored in `genetox_smarts`
        if any(mol.HasSubstructMatch(e) for v in genetox_smarts.values() for e in v):
            return 'Genetox'

        return False  # No genetox alerts found

    def decision_tree(self):
    # Check if the substance is inorganic
        inorg_filter = (
            self.metal_ions() != 'essential metal ion' and 
            (self.has_metal_atom() == 'inorganic' or self.P_inorg() == 'inorganic')
        )

        if inorg_filter:
            return 'Inorganic - TTC not applicable'

        # Check for DLC, Steroid, or HPC
        if self.dlc() == 'Dioxin-like':
            return 'DLC - TTC not applicable'
        if self.steroid() == 'Steroid':
            return 'Steroid - TTC not applicable'
        if self.hpc() == 'HPC':
            return 'HPC - TTC not applicable'

        # Check for genetic toxicity *before* OPCs
        if self.genetox_alerts() == 'Genetox':
            return 'GeneTox - TTC of 0.15 ug/day'

        # Check for OPCs
        if self.opc() == 'OP or carbamate':
            return 'OP_carbamate - TTC of 18 ug/day'

        # If none of the above conditions were met, apply Cramer classification
        return 'Likely Cramer class applicable'    
        
        
   
    


    def process_batch(self, df):
        results = []
        for index, row in df.iterrows():
            try:
                self.add_chemical(row['dtxsid'], row['smiles'])
                result = self.decision_tree()
            except ValueError as e:
                result = str(e)
            results.append(result)
        df['final_TTC_category'] = results
        return df
        




 