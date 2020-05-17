from calculations.python.paths import *
import csv

def get_genes() -> [str]:
    genes = []
    with open(data_path(REFERENCE_GENES), 'r') as file:
        for line in file:
            genes.append(line.split()[0])
    return genes

genes = get_genes()
comparison_columns = ["compared", "distance"]
genome_columns = ["header"]
gene_columns = ["begin", "end", "invalid_nucleotides"]
gene_comparison_columns = ["mutations", "protein_mutations"]

def entry_to_csv_dict(entry):
    result = {}
    for col in comparison_columns:
        result[col] = entry[1][col]
    for col in genome_columns:
        result[col] = entry[0][col]
    result["date"] = result["compared"].split('|')[2]
    if result["distance"] < 0:
        return None
    genes_dict = {g["name"]: g for g in entry[0]["genes"]}
    genes_comparison_dict = {g["name"]: g for g in entry[1]["genes"]}
    for gene in genes:
        if gene not in genes_dict or gene not in genes_comparison_dict:
            return None
        for col in gene_columns:
            result[f"{gene}_{col}"] = genes_dict[gene][col]
        for col in gene_comparison_columns:
            result[f"{gene}_{col}"] = [(mutation["position"], '?', mutation["arg"]) for mutation in genes_comparison_dict[gene][col]]
        result[f"{gene}_distance"] = len(genes_comparison_dict[gene]["mutations"])
        result[f"{gene}_protein_distance"] = len(genes_comparison_dict[gene]["protein_mutations"])
    return result

def get_csv_columns() -> [str]:
    return [f"{gene}_{col}" for gene in genes for col in gene_columns] + [f"{gene}_{col}" for gene in genes for col in ["distance", "protein_distance"] + gene_comparison_columns] + comparison_columns + ["date"] + genome_columns

with open(data_path(CALCULATED_GENOMES_DATA), 'r') as input, open(data_path(DISTANCES_CSV), 'w') as output:
    next(input)
    writer = csv.DictWriter(output, fieldnames=get_csv_columns())
    writer.writeheader()
    errors = 0
    total = 0
    for line in input:
        row = entry_to_csv_dict(eval(line))
        total += 1
        if row is not None:
            writer.writerow(row)
        else:
            errors += 1
    print(f"{errors} erroneous entries out of {total}")