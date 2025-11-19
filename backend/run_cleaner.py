import argparse
import csv
from src.cleaner import clean_molecule

def process_file(input_path, output_path):
    total = 0
    cleaned_count = 0

    with open(input_path, 'r') as infile, open(output_path, 'w', newline='') as outfile:
        reader = (line.strip() for line in infile if line.strip())
        writer = csv.writer(outfile)
        writer.writerow(['original_smiles', 'cleaned_smiles'])

        for line in reader:
            total += 1
            # If tab/space separated name, split
            parts = line.split()
            original_smiles = parts[0]
            cleaned = clean_molecule(original_smiles)

            if cleaned:
                writer.writerow([original_smiles, cleaned])
                cleaned_count += 1
            else:
                writer.writerow([original_smiles, 'CLEANING_FAILED'])

            if total % 100 == 0:
                print(f"Processed {total} molecules...")

    print(f"Finished processing {total} molecules.")
    print(f"Successfully cleaned: {cleaned_count}")
    print(f"Failed: {total - cleaned_count}")

def main():
    parser = argparse.ArgumentParser(description="Clean SMILES file using RDKit pipeline")
    parser.add_argument('--input', '-i', required=True, help="Input SMILES file path (one per line)")
    parser.add_argument('--output', '-o', required=True, help="Output CSV file path")

    args = parser.parse_args()
    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
