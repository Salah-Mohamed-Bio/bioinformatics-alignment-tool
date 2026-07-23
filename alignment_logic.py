def global_alignment(seq1, seq2, MATCH, MISMATCH, GAP):
    
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    Matrix = [[0 for j in range(cols)] for i in range(rows)]
    # traceback pointer matrix: 'd' diagonal, 'u' up, 'l' left
    Pointer = [[None for j in range(cols)] for i in range(rows)]

    # first column
    for i in range(rows):
        Matrix[i][0] = i * GAP
        if i > 0:
            Pointer[i][0] = 'u'

    # first row
    for j in range(cols):
        Matrix[0][j] = j * GAP
        if j > 0:
            Pointer[0][j] = 'l'


    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = Matrix[i-1][j-1] + (MATCH if seq1[i-1] == seq2[j-1] else MISMATCH)
            up = Matrix[i-1][j] + GAP
            left = Matrix[i][j-1] + GAP

            best = max(diagonal, up, left)
            Matrix[i][j] = best

            if best == diagonal:
                Pointer[i][j] = 'd'
            elif best == up:
                Pointer[i][j] = 'u'
            else:
                Pointer[i][j] = 'l'


    align1 = ""
    align2 = ""

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and Pointer[i][j] == 'd':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
            
        elif i > 0 and Pointer[i][j] == 'u':
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:  # Pointer == 'l' (or j only)
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1


    return align1, align2,  Matrix [len(seq1)][len(seq2)] 


def local_alignment(seq1, seq2, MATCH, MISMATCH, GAP):
    
    rows = len(seq1) + 1
    cols = len(seq2) + 1


    Matrix = [[0 for j in range(cols)] for i in range(rows)]
    # traceback pointer matrix: 'd' diagonal, 'u' up, 'l' left, '0' stop
    Pointer = [[None for j in range(cols)] for i in range(rows)]


    max_score = 0
    max_position = (0, 0)


    for i in range(1, rows):
        for j in range(1, cols):

            score = MATCH if seq1[i-1] == seq2[j-1] else MISMATCH

            diagonal = Matrix[i-1][j-1] + score
            up = Matrix[i-1][j] + GAP
            left = Matrix[i][j-1] + GAP

            best = max(0, diagonal, up, left)
            Matrix[i][j] = best

            # store the source that produced the value (prefer 0/stop when best == 0)
            if best == 0:
                Pointer[i][j] = '0'
            elif best == diagonal:
                Pointer[i][j] = 'd'
            elif best == up:
                Pointer[i][j] = 'u'
            else:
                Pointer[i][j] = 'l'

            if Matrix[i][j] > max_score:
                max_score = Matrix[i][j]
                max_position = (i, j)


    align1 = ""
    align2 = ""

    # start backtrack from the max-scoring cell, not bottom-right
    i, j = max_position

    while Pointer[i][j] != '0' and i > 0 and j > 0:

        if Pointer[i][j] == 'd':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif Pointer[i][j] == 'u':
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:  # 'l'
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    
    return align1, align2, max_score 



def get_symbol(align1, align2):

    symbol = ""
    for a, b in zip(align1, align2):
        if a == b:
            symbol += "|"
        elif a == "-" or b == "-":
            symbol += " "
        else:
            symbol += "."
    return symbol


# ---------------------------------------------------------------------------
# Neighbor-Joining (NJ) Phylogenetic Tree Construction
# ---------------------------------------------------------------------------

def neighbor_joining(distance_matrix, labels):
    """
    Construct a phylogenetic tree using the Neighbor-Joining algorithm.
    
    Parameters:
        distance_matrix: list of lists - symmetric distance matrix (n x n)
        labels: list of strings - sequence/OTU labels
        
    Returns:
        Newick format tree string
    """
    import copy
    n = len(labels)
    if n < 2:
        return ""
    if n == 2:
        return f"({labels[0]}:{distance_matrix[0][1]/2:.4f},{labels[1]}:{distance_matrix[0][1]/2:.4f})"
    
    # Convert to float and work with copies
    # Pre-allocate a larger matrix that can hold internal nodes (max size: n + n - 2 = 2n)
    max_size = 2 * n
    D = [[0.0] * max_size for _ in range(max_size)]
    for i in range(n):
        for j in range(n):
            D[i][j] = distance_matrix[i][j]
    
    # Work with a list of current (remaining) taxa indices
    taxa = list(range(n))
    # Track the names of active nodes
    active = {i: labels[i] for i in range(n)}
    
    # Newick tree pieces: for each node, store (left_child, right_child, left_dist, right_dist)
    tree = {}
    next_id = n  # For naming internal nodes
    
    while len(taxa) > 2:
        m = len(taxa)
        
        # Step 1: Calculate total distances
        total_dist = {}
        for i in taxa:
            total_dist[i] = sum(D[i][j] for j in taxa)
        
        # Step 2: Calculate Q matrix (for this subset)
        min_q = float('inf')
        min_pair = None
        
        for idx_i, i in enumerate(taxa):
            for idx_j, j in enumerate(taxa):
                if i >= j:
                    continue
                q_val = (m - 2) * D[i][j] - total_dist[i] - total_dist[j]
                if q_val < min_q:
                    min_q = q_val
                    min_pair = (i, j)
        
        i, j = min_pair
        
        # Step 3: Calculate distances from new node u to i and j
        d_iu = 0.5 * D[i][j] + (total_dist[i] - total_dist[j]) / (2 * (m - 2))
        d_ju = D[i][j] - d_iu
        
        # Step 4: Create new node u
        u = next_id
        next_id += 1
        
        # Store tree info
        tree[u] = {
            'left': i,
            'right': j,
            'left_dist': round(d_iu, 4),
            'right_dist': round(d_ju, 4),
            'name': f"Internal{u}"
        }
        
        # Step 5: Calculate distances from u to all other active nodes
        for k in taxa:
            if k != i and k != j:
                D[u][k] = 0.5 * (D[i][k] + D[j][k] - D[i][j])
                D[k][u] = D[u][k]
        
        # Step 6: Replace i and j with u
        taxa.remove(j)
        taxa.remove(i)
        taxa.append(u)
        active[u] = f"Internal{u}"
    
    # Final step: connect the last two nodes
    if len(taxa) == 2:
        i, j = taxa
        dist = D[i][j]
        # Return Newick string
        newick = _build_newick(i, j, dist, tree, active)
        return newick
    
    return ""


def _build_newick(i, j, dist, tree, active):
    """Build Newick string recursively."""
    def get_subtree(node, branch_length=0):
        # Check internal nodes first
        if node in tree:
            left = tree[node]['left']
            right = tree[node]['right']
            left_dist = tree[node]['left_dist']
            right_dist = tree[node]['right_dist']
            left_str = get_subtree(left, left_dist)
            right_str = get_subtree(right, right_dist)
            return f"({left_str},{right_str}):{branch_length:.4f}"
        
        # Leaf node
        name = active.get(node, f"Taxon{node}")
        return f"{name}:{branch_length:.4f}"
    
    # Handle the final two nodes
    left_str = get_subtree(i, dist/2)
    right_str = get_subtree(j, dist/2)
    return f"({left_str},{right_str});"


def convert_score_to_distance(similarity_matrix):
    """
    Convert a similarity score matrix to a distance matrix.
    Uses: distance = (max_score - score) / max_score  (normalized)
    """
    n = len(similarity_matrix)
    if n == 0:
        return []
    
    # Find max and min
    max_val = max(max(row) for row in similarity_matrix)
    min_val = min(min(row) for row in similarity_matrix)
    range_val = max_val - min_val if max_val != min_val else 1
    
    dist = []
    for i in range(n):
        row = []
        for j in range(n):
            # Normalize distance: 1 - (score - min) / range
            d = 1.0 - (similarity_matrix[i][j] - min_val) / range_val
            if i == j:
                d = 0.0  # same sequence = zero distance
            row.append(round(d, 6))
        dist.append(row)
    return dist


if __name__=="__main__":

    seq1 = input("enter first sequenec: ").strip().upper()
    seq2 = input("enter the second sequence: ").strip().upper()
    
    MATCH = int(input("enter Match SCORE: "))
    MISMATCH = int(input("enter Mismatch SCORE: "))
    GAP = int(input("enter THE Gap: "))

    choice = input("choice 1 or 2 ; 1:Global, 2:Local:")

    if choice == '1':
        align1, align2, score = global_alignment(seq1, seq2, MATCH, MISMATCH, GAP)

    else:
        align1, align2, score = local_alignment(seq1, seq2, MATCH, MISMATCH, GAP)



    def print_alignment(align1, align2, symbol, width=100):
        for i in range(0, len(align1), width):
            print(align1[i:i+width])
            print(symbol[i:i+width])
            print(align2[i:i+width])
            print()
    symbol = get_symbol(align1, align2)    
    print("\nAlignment")
    print_alignment(align1, align2, symbol)
    print("\nscore =",score)