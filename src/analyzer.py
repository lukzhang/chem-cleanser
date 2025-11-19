from rdkit import Chem
from rdkit.Chem import Descriptors

def calc_descriptors(mol):
    return {
        'MW': round(Descriptors.MolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'TPSA': round(Descriptors.TPSA(mol), 2),
        'RotatableBonds': Descriptors.NumRotatableBonds(mol),
    }

def lipinski_rule_of_5(desc):
    return (desc['MW'] <= 500 and
            desc['LogP'] <= 5 and
            desc['HBD'] <= 5 and
            desc['HBA'] <= 10)

def veber_rule(desc):
    return (desc['RotatableBonds'] <= 10 and desc['TPSA'] <= 140)
