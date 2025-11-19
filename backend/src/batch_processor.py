def read_smiles_file(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                yield line.split()[0]  # return just the SMILES part

def write_results_csv(filepath, rows, header):
    import csv
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
