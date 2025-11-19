from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

remover = SaltRemover.SaltRemover()
normalizer = rdMolStandardize.Normalizer()
reionizer = rdMolStandardize.Reionizer()
tautomer_enumerator = rdMolStandardize.TautomerEnumerator()

def clean_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = remover.StripMol(mol, dontRemoveEverything=True)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    mol = normalizer.normalize(mol)
    mol = reionizer.reionize(mol)
    mol = tautomer_enumerator.Canonicalize(mol)

    # Perform sanitization
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)

    if mol is None or mol.GetNumAtoms() == 0:
        return None

    can_smiles = Chem.MolToSmiles(mol, canonical=True)
    return can_smiles
