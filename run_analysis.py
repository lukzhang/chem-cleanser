import argparse
import csv
from src.cleaner import clean_molecule
from src.analyzer import calc_descriptors, lipinski_rule_of_5, veber_rule
from rdkit import Chem

def process_file(input_path, output_path):
    total = 0
    cleaned_count = 0

    with open(input_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
        reader = (line.strip() for line in infile if line.strip())
        writer = csv.writer(outfile)

        # Header with descriptors and rules
        writer.writerow([
            'original_smiles', 'cleaned_smiles', 'MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'RotatableBonds',
            'Lipinski_pass', 'Veber_pass'
        ])

        for line in reader:
            total += 1
            parts = line.split()
            original_smiles = parts[0]
            cleaned_smiles = clean_molecule(original_smiles)

            if cleaned_smiles:
                mol = Chem.MolFromSmiles(cleaned_smiles)
                if mol:
                    desc = calc_descriptors(mol)
                    lipinski = lipinski_rule_of_5(desc)
                    veber = veber_rule(desc)

                    writer.writerow([
                        original_smiles, cleaned_smiles,
                        desc['MW'], desc['LogP'], desc['HBD'], desc['HBA'], desc['TPSA'], desc['RotatableBonds'],
                        lipinski, veber
                    ])
                    cleaned_count += 1
                else:
                    writer.writerow([original_smiles, 'INVALID_AFTER_CLEANING'] + [''] * 8)
            else:
                writer.writerow([original_smiles, 'CLEANING_FAILED'] + [''] * 8)

            if total % 100 == 0:
                print(f"Processed {total} molecules...")

    print(f"Finished processing {total} molecules.")
    print(f"Successfully analyzed: {cleaned_count}")
    print(f"Failed: {total - cleaned_count}")

def main():
    parser = argparse.ArgumentParser(description="Clean and analyze SMILES for drug-likeness")
    parser.add_argument('--input', '-i', required=True, help="Input SMILES file path")
    parser.add_argument('--output', '-o', required=True, help="Output CSV file path")

    args = parser.parse_args()
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
