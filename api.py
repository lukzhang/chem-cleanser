# api.py

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List
from rdkit import Chem

from src.cleaner import clean_molecule
from src.analyzer import calc_descriptors, lipinski_rule_of_5, veber_rule

from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="Chem-Cleanser API", version="0.1.0")

class SmilesList(BaseModel):
    smiles: List[str]

class AnalyzeResponse(BaseModel):
    smiles: List[str]
    lipinski_pass: List[bool]
    veber_pass: List[bool]

@app.post("/clean")
def clean(smiles_data: SmilesList):
    cleaned = []
    for smi in smiles_data.smiles:
        cleaned_smi = clean_molecule(smi)
        cleaned.append(cleaned_smi)
    return {"cleaned_smiles": cleaned}

@app.post("/analyze", response_model=AnalyzeResponse)
def analyze_smiles(request: SmilesList):
    results, lipinski, veber = analyze_molecules(request.smiles)
    return AnalyzeResponse(
        smiles=results,
        lipinski_pass=lipinski,
        veber_pass=veber
    )

@app.post("/pipeline", response_model=AnalyzeResponse)
def pipeline(smiles_data: SmilesList):
    try:
        # clean + analyze combined
        cleaned = [clean_molecule(s) for s in smiles_data.smiles]
        results, lipinski, veber = analyze_molecules(cleaned)
        return AnalyzeResponse(
            smiles=results,
            lipinski_pass=lipinski,
            veber_pass=veber
        )
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

def analyze_molecules(smiles_list: List[str]):
    results = []
    lipinski = []
    veber = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append("")  # or "Invalid SMILES"
            lipinski.append(False)
            veber.append(False)
            continue
        desc = calc_descriptors(mol)
        results.append(smi)
        lipinski.append(lipinski_rule_of_5(desc))
        veber.append(veber_rule(desc))
    return results, lipinski, veber


# CORS setup
origins = [
    "http://localhost",
    "http://localhost:8000",
    "http://127.0.0.1",
    "http://127.0.0.1:8000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
