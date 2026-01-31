def process_bipartite(input_path, output_path):
    left_map = {}
    right_map = {}
    edges_set = set()

    with open(input_path, "r") as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("%")]

    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        u, v = int(parts[0]), int(parts[1])

        # rinumera lato sinistro
        if u not in left_map:
            left_map[u] = len(left_map)

        # rinumera lato destro
        if v not in right_map:
            right_map[v] = len(right_map)

        # evita duplicati
        edges_set.add((left_map[u], right_map[v]))

    n_left = len(left_map)
    n_right = len(right_map)
    n_edges = len(edges_set)

    with open(output_path, "w") as f:
        f.write(f"{n_left}\n")
        f.write(f"{n_right}\n")
        f.write(f"{n_edges}\n")
        for u, v in edges_set:
            # lato destro parte da n_left
            f.write(f"{u} {n_left + v}\n")

if __name__ == "__main__":
    input_file = "./dataset/raw/power.tsv"
    output_file = "./dataset/data/power.edges"
    process_bipartite(input_file, output_file)
