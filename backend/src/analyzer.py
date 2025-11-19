from rdkit.Chem import Descriptors

def calc_descriptors(mol):
    desc = {}
    desc['MW'] = Descriptors.MolWt(mol)
    desc['LogP'] = Descriptors.MolLogP(mol)
    desc['HBD'] = Descriptors.NumHDonors(mol)
    desc['HBA'] = Descriptors.NumHAcceptors(mol)
    desc['TPSA'] = Descriptors.TPSA(mol)
    desc['RotatableBonds'] = Descriptors.NumRotatableBonds(mol)
    return desc

def lipinski_rule_of_5(desc):
    if desc['MW'] > 500:
        return False
    if desc['LogP'] > 5:
        return False
    if desc['HBD'] > 5:
        return False
    if desc['HBA'] > 10:
        return False
    return True

def veber_rule(desc):
    if desc['TPSA'] > 140:
        return False
    if desc['RotatableBonds'] > 10:
        return False
    return True

def molecular_weight_rule(desc):
    return 150 <= desc['MW'] <= 500
