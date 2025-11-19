# api.py

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, validator, Field
from typing import List
from rdkit import Chem

from src.cleaner import clean_molecule
from src.analyzer import calc_descriptors, lipinski_rule_of_5, veber_rule

from src.analyzer import calc_descriptors, lipinski_rule_of_5, veber_rule, molecular_weight_rule

from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="Chem-Cleanser API", version="0.1.0")

class SmilesList(BaseModel):
    smiles: List[str] = Field(
        ...,
        example=["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    )

    @validator('smiles')
    def check_not_empty(cls, v):
        if not v or len(v) == 0:
            raise ValueError('smiles list must not be empty')
        return v

class AnalyzeResponse(BaseModel):
    smiles: List[str]
    lipinski_pass: List[bool]
    veber_pass: List[bool]
    mw_pass: List[bool]  # Molecular Weight filter pass/fail


@app.post("/clean")
def clean(smiles_data: SmilesList):
    cleaned = []
    for smi in smiles_data.smiles:
        cleaned_smi = clean_molecule(smi)
        cleaned.append(cleaned_smi)
    return {"cleaned_smiles": cleaned}

@app.post("/analyze", response_model=AnalyzeResponse)
def analyze_smiles(request: SmilesList):
    results, lipinski, veber, mw = analyze_molecules(request.smiles)
    return AnalyzeResponse(
        smiles=results,
        lipinski_pass=lipinski,
        veber_pass=veber,
        mw_pass=mw
    )

@app.post("/pipeline", response_model=AnalyzeResponse)
def pipeline(smiles_data: SmilesList):
    try:
        cleaned = [clean_molecule(s) for s in smiles_data.smiles]
        results, lipinski, veber, mw = analyze_molecules(cleaned)
        return AnalyzeResponse(
            smiles=results,
            lipinski_pass=lipinski,
            veber_pass=veber,
            mw_pass=mw
        )
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

def analyze_molecules(smiles_list: List[str]):
    results = []
    lipinski = []
    veber = []
    mw = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise HTTPException(status_code=400, detail=f"Invalid SMILES string: {smi}")
        desc = calc_descriptors(mol)
        results.append(smi)
        lipinski.append(lipinski_rule_of_5(desc))
        veber.append(veber_rule(desc))
        mw.append(molecular_weight_rule(desc))
    return results, lipinski, veber, mw


# CORS setup
origins = [
    "http://localhost",
    "http://localhost:8000",
    "http://127.0.0.1",
    "http://127.0.0.1:8000",
    "http://localhost:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
