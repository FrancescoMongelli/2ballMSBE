def convert_drug_gene(input_path, output_path):
    drug_map = {}
    gene_map = {}
    edges = set()

    with open(input_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            drug, gene = line.split()

            if drug not in drug_map:
                drug_map[drug] = len(drug_map)

            if gene not in gene_map:
                gene_map[gene] = len(gene_map)

            edges.add((drug_map[drug], gene_map[gene]))

    n_left = len(drug_map)
    n_right = len(gene_map)
    n_edges = len(edges)

    with open(output_path, "w") as f:
        f.write(f"{n_left}\n")
        f.write(f"{n_right}\n")
        f.write(f"{n_edges}\n")

        for u, v in edges:
            f.write(f"{u} {n_left + v}\n")


if __name__ == "__main__":
    convert_drug_gene(
        "./dataset/raw/ChG.tsv",
        "./dataset/data/ChG.edges"
    )
