seq1 = input("enter the first sequence: ")
seq2 = input("enter the second sequence: ")

if not seq1 or not seq2:
    raise SystemExit("Error: both sequences must be non-empty.")

MATCH = int(input("enter MATCH SCORE: "))
MISMATCH = int(input("enter MISMATCH SCORE: "))
GAP = int(input("enter THE GAP: "))


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

while Pointer[i][j] != '0':
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


symbol = ""
for a, b in zip(align1, align2):
    if a == b:
        symbol += "|"
    elif a == "-" or b == "-":
        symbol += " "
    else:
        symbol += "."


def print_alignment(align1, align2, symbol, width=100):
    for i in range(0, len(align1), width):
        print(align1[i:i+width])
        print(symbol[i:i+width])
        print(align2[i:i+width])
        print()


print("\nAlignment")
print_alignment(align1, align2, symbol)
print("\nmax_score =", max_score)
