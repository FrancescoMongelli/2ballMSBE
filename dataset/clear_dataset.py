left_nodes = set()
right_nodes = set()
edges = set()

with open("./dataset/raw/it.csv","r") as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue

        src, trg, _, _ = line.strip().split(",")
        src = int(src)
        trg = int(trg)

        left_nodes.add(src)
        right_nodes.add(trg)
        edges.add((src, trg))

# tutti i nodi dichiarati
all_declared_nodes = set(range(max(right_nodes)+1))

# nodi che non compaiono negli archi
isolated = all_declared_nodes - (left_nodes.union(right_nodes))

print("Nodi isolati:", isolated)

while isolated:
    iso = isolated.pop()
    # esempio di regola semplice:
    # se l'ID isolato <= max(left_nodes) -> è un left
    if iso <= max(left_nodes):
        left_nodes.add(iso)
    else:
        right_nodes.add(iso)

num_left = len(left_nodes)
num_right = len(right_nodes)
num_edges = len(edges)

with open("./dataset/data/it.edges", "w") as out:
    out.write(f"{num_left}\n")
    out.write(f"{num_right}\n")
    out.write(f"{num_edges}\n")
    
    for src, trg in sorted(edges):
        out.write(f"{src} {trg}\n")
