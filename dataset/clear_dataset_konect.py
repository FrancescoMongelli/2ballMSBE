def convert_file(input_path, output_path):
    with open(input_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Salta le prime tre righe (header)
    edges = lines[3:]

    left_map = {}
    right_map = {}
    next_left = 0
    next_right = 0

    # Prima passata: mappiamo tutti i nodi per ottenere il totale
    for line in edges:
        l, r = map(int, line.split())
        if l not in left_map:
            left_map[l] = next_left
            next_left += 1
        if r not in right_map:
            right_map[r] = next_right
            next_right += 1

    n_left = next_left
    n_right = next_right

    # Seconda passata: scriviamo gli archi con shift dei nodi di destra
    new_edges = []
    for line in edges:
        l, r = map(int, line.split())
        new_edges.append(f"{left_map[l]} {n_left + right_map[r]}")  # shift right nodes

    # Scrittura output: prime tre righe = n_left, n_right, num_arcs
    with open(output_path, "w") as f:
        f.write(f"{n_left}\n")
        f.write(f"{n_right}\n")
        f.write(f"{len(new_edges)}\n")
        for line in new_edges:
            f.write(line + "\n")


if __name__ == "__main__":
    convert_file("./dataset/raw/AC", "./dataset/data/AC.edges")
