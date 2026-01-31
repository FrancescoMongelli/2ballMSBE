def remap_edges(input_path, output_path):
    left_map = {}   # mapping dei nodi di sinistra
    right_map = {}  # mapping dei nodi di destra
    edges = set()   # insieme per rimuovere duplicati

    with open(input_path, "r") as f:
        for line in f:
            if line.startswith("%") or line.strip() == "":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            u, v = int(parts[0]), int(parts[1])

            if u not in left_map:
                left_map[u] = len(left_map)
            if v not in right_map:
                right_map[v] = len(right_map)

            u_new = left_map[u]
            v_new = right_map[v]

            # usa una tupla per rimuovere duplicati
            edges.add((u_new, v_new))

    n_left = len(left_map)
    n_right = len(right_map)
    n_edges = len(edges)

    with open(output_path, "w") as out:
        out.write(f"{n_left}\n")
        out.write(f"{n_right}\n")
        out.write(f"{n_edges}\n")
        for u_new, v_new in edges:
            out.write(f"{u_new} {v_new + n_left}\n")  # shift nodi di destra

    print(f"Done. Left nodes: {n_left}, Right nodes: {n_right}, Edges: {n_edges}")

# Esempio d'uso
# remap_edges("input_file.txt", "output_file.txt")
if __name__ == "__main__":
    remap_edges("./dataset/raw/AM", "./dataset/data/AM.edges")