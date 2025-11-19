FROM continuumio/miniconda3

WORKDIR /app

# Install RDKit from conda-forge
RUN conda install -c conda-forge rdkit -y

# Copy your code
COPY ./src ./src
COPY run_analysis.py .
COPY run_cleaner.py .

CMD ["python", "run_analysis.py", "--input", "data/sample_smiles.txt", "--output", "data/analysis_output.csv"]
