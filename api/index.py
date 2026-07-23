"""
Vercel Python Serverless Function for Sequence Alignment.
"""
import sys
import os
import io
import re
from flask import Flask, request, jsonify, render_template

# Add parent directory to path so we can import alignment_logic
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from alignment_logic import global_alignment, local_alignment, get_symbol, neighbor_joining, convert_score_to_distance

app = Flask(__name__)

# Tell Flask where to find templates
app.template_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

# Vercel free tier has ~4.5MB request limit
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024
# Maximum sequence length (increase as needed, but watch memory)
MAX_SEQ_LENGTH = 10000


# Manual CORS
@app.after_request
def add_cors_headers(response):
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, OPTIONS"
    return response


def parse_fasta_text(text):
    """Parse a FASTA string (without BioPython) and return the sequence."""
    lines = text.strip().splitlines()
    seq_parts = []
    for line in lines:
        if line.startswith(">"):
            continue
        clean = re.sub(r"[^A-Za-z]", "", line)
        if clean:
            seq_parts.append(clean)
    if not seq_parts:
        raise ValueError("No sequence data found in FASTA text.")
    return "".join(seq_parts).upper()


def parse_fasta_file(file_storage):
    """Parse a FASTA file uploaded via multipart form."""
    try:
        text = file_storage.read().decode("utf-8")
        return parse_fasta_text(text)
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA file: {e}")


def parse_multi_fasta(text):
    """Parse a multi-FASTA string into a list of (id, sequence) tuples."""
    records = []
    current_id = None
    current_seq_parts = []
    for line in text.strip().splitlines():
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, "".join(current_seq_parts).upper()))
            current_id = line[1:].strip().split()[0]
            current_seq_parts = []
        else:
            clean = re.sub(r"[^A-Za-z]", "", line)
            if clean:
                current_seq_parts.append(clean)
    if current_id is not None:
        records.append((current_id, "".join(current_seq_parts).upper()))
    return records


def run_alignment(seq1, seq2, match, mismatch, gap, mode):
    """Run global or local alignment and return NCBI-style results."""
    if not seq1 or not seq2:
        raise ValueError("Both sequences are required.")
    if len(seq1) > MAX_SEQ_LENGTH or len(seq2) > MAX_SEQ_LENGTH:
        raise ValueError(
            f"Sequences are too long (max {MAX_SEQ_LENGTH} bp). "
            f"Got {len(seq1)} and {len(seq2)} bp."
        )
    if mode == "global":
        align1, align2, score = global_alignment(seq1, seq2, match, mismatch, gap)
    elif mode == "local":
        align1, align2, score = local_alignment(seq1, seq2, match, mismatch, gap)
    else:
        raise ValueError("mode must be 'global' or 'local'.")
    symbol = get_symbol(align1, align2)
    
    # Compute NCBI-style statistics
    aligned_len = len(align1)
    identities = sum(1 for a, b in zip(align1, align2) if a == b)
    gaps = sum(1 for a, b in zip(align1, align2) if a == '-' or b == '-')
    mismatches = aligned_len - identities - gaps
    identity_pct = round(100 * identities / aligned_len, 1) if aligned_len > 0 else 0
    gap_pct = round(100 * gaps / aligned_len, 1) if aligned_len > 0 else 0
    mismatch_pct = round(100 * mismatches / aligned_len, 1) if aligned_len > 0 else 0
    
    return {
        "align1": align1,
        "align2": align2,
        "symbol": symbol,
        "score": score,
        "mode": mode,
        "seq1_original": seq1,
        "seq2_original": seq2,
        "aligned_length": aligned_len,
        "identities": identities,
        "identity_pct": identity_pct,
        "mismatches": mismatches,
        "mismatch_pct": mismatch_pct,
        "gaps": gaps,
        "gap_pct": gap_pct,
    }


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------


@app.route("/", methods=["GET"])
def index():
    """Serve the web frontend."""
    return render_template("index.html")


@app.route("/health", methods=["GET"])
def health():
    """Health check endpoint."""
    return jsonify({"status": "ok"})


@app.route("/align", methods=["POST"])
def align():
    """
    Run sequence alignment from JSON payload.
    """
    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Request body must be valid JSON."}), 400

    seq1 = data.get("seq1", "").strip().upper()
    seq2 = data.get("seq2", "").strip().upper()
    match = data.get("match", 2)
    mismatch = data.get("mismatch", -3)
    gap = data.get("gap", -5)
    mode = data.get("mode", "global")

    try:
        result = run_alignment(seq1, seq2, match, mismatch, gap, mode)
        return jsonify(result), 200
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except MemoryError:
        return jsonify({"error": "Not enough memory for this alignment. Please use shorter sequences."}), 400


@app.route("/align/fasta", methods=["POST"])
def align_fasta():
    """
    Run sequence alignment from uploaded FASTA files.
    """
    if "file1" not in request.files or "file2" not in request.files:
        return jsonify({"error": "Both 'file1' and 'file2' are required."}), 400

    try:
        seq1 = parse_fasta_file(request.files["file1"])
        seq2 = parse_fasta_file(request.files["file2"])
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    match = request.form.get("match", 2, type=int)
    mismatch = request.form.get("mismatch", -3, type=int)
    gap = request.form.get("gap", -5, type=int)
    mode = request.form.get("mode", "global")

    try:
        result = run_alignment(seq1, seq2, match, mismatch, gap, mode)
        return jsonify(result), 200
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except MemoryError:
        return jsonify({"error": "Not enough memory for this alignment. Please use shorter sequences."}), 400


@app.route("/align/batch", methods=["POST"])
def align_batch():
    """
    Run batch all-vs-all alignment from a multi-FASTA file upload.
    """
    import time
    start_time = time.time()

    if "file" not in request.files:
        return jsonify({"error": "Multi-FASTA file is required."}), 400

    try:
        text = request.files["file"].read().decode("utf-8")
        records = parse_multi_fasta(text)
    except Exception as e:
        return jsonify({"error": f"Failed to parse file: {e}"}), 400

    if len(records) < 2:
        return jsonify({"error": "Multi-FASTA file must contain at least 2 sequences."}), 400

    if len(records) > 20:
        return jsonify({"error": "Maximum 20 sequences allowed for batch alignment."}), 400

    match = request.form.get("match", 2, type=int)
    mismatch = request.form.get("mismatch", -3, type=int)
    gap = request.form.get("gap", -5, type=int)

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    n = len(seqs)

    # Check if sequences are too long for batch
    for i, seq in enumerate(seqs):
        if len(seq) > 3000:
            return jsonify({
                "error": f"Sequence '{ids[i]}' is {len(seq)} bp long. "
                         f"Vercel free tier limit: 3000 bp per sequence for batch alignment. "
                         f"Use pairwise mode or run locally."
            }), 400

    # Estimate if it will timeout (n sequences × n comparisons × seq_length²)
    # Rough estimate: each comparison of ~1000bp takes ~0.1s on Vercel
    estimated_pairs = (n * (n + 1)) / 2
    estimated_time = estimated_pairs * (max(len(s) for s in seqs) / 1000) ** 2 * 0.03
    if estimated_time > 8:
        return jsonify({
            "error": f"Too many/long sequences for Vercel's 10-second timeout. "
                     f"Estimated time: {estimated_time:.1f}s. "
                     f"Please use fewer sequences (max {n-1}) or shorter ones, or run locally."
        }), 400

    matrix = []
    try:
        for i in range(n):
            row = []
            for j in range(n):
                if j < i:
                    row.append(matrix[j][i])
                else:
                    # Check remaining time
                    elapsed = time.time() - start_time
                    if elapsed > 8:
                        raise TimeoutError("Vercel 10-second timeout approaching. Please use fewer/shorter sequences.")
                    _, _, score = global_alignment(seqs[i], seqs[j], match, mismatch, gap)
                    row.append(score)
            matrix.append(row)
    except MemoryError:
        return jsonify({"error": "Not enough memory. Please use shorter sequences or fewer sequences."}), 400
    except TimeoutError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": f"Alignment error: {str(e)[:200]}"}), 400

    return jsonify({
        "ids": ids,
        "matrix": matrix,
    }), 200


@app.route("/phylogeny", methods=["POST"])
def phylogeny():
    """
    Generate a phylogenetic tree from a multi-FASTA file using Neighbor-Joining.
    
    Accepts multipart form with:
        file - multi-FASTA file (2-15 sequences)
        match - match score (default: 5)
        mismatch - mismatch score (default: -1)
        gap - gap penalty (default: -3)
    
    Returns:
        newick - Newick format tree string
        ids - sequence labels
        distance_matrix - normalized distance matrix
        similarity_matrix - raw similarity matrix
    """
    import time
    start_time = time.time()

    if "file" not in request.files:
        return jsonify({"error": "Multi-FASTA file is required."}), 400

    try:
        text = request.files["file"].read().decode("utf-8")
        records = parse_multi_fasta(text)
    except Exception as e:
        return jsonify({"error": f"Failed to parse file: {e}"}), 400

    if len(records) < 3:
        return jsonify({"error": "Phylogenetic tree requires at least 3 sequences."}), 400

    if len(records) > 15:
        return jsonify({"error": "Maximum 15 sequences allowed for tree building."}), 400

    match = request.form.get("match", 2, type=int)
    mismatch = request.form.get("mismatch", -3, type=int)
    gap = request.form.get("gap", -5, type=int)

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    n = len(seqs)

    # Check sequence lengths
    for i, seq in enumerate(seqs):
        if len(seq) > 2000:
            return jsonify({
                "error": f"Sequence '{ids[i]}' is {len(seq)} bp long. "
                         f"Maximum 2000 bp for tree building."
            }), 400

    # Build similarity matrix
    similarity_matrix = []
    try:
        for i in range(n):
            row = []
            for j in range(n):
                elapsed = time.time() - start_time
                if elapsed > 8:
                    raise TimeoutError("Vercel timeout approaching. Please use fewer/shorter sequences.")
                if j < i:
                    row.append(similarity_matrix[j][i])
                else:
                    _, _, score = global_alignment(seqs[i], seqs[j], match, mismatch, gap)
                    row.append(score)
            similarity_matrix.append(row)
    except TimeoutError as e:
        return jsonify({"error": str(e)}), 400
    except MemoryError:
        return jsonify({"error": "Not enough memory for this alignment."}), 400
    except Exception as e:
        return jsonify({"error": f"Alignment error: {str(e)[:200]}"}), 400

    # Convert to distance matrix
    distance_matrix = convert_score_to_distance(similarity_matrix)

    # Build neighbor-joining tree
    try:
        newick = neighbor_joining(distance_matrix, ids)
    except Exception as e:
        return jsonify({"error": f"Tree building error: {str(e)[:200]}"}), 400

    return jsonify({
        "newick": newick,
        "ids": ids,
        "distance_matrix": distance_matrix,
        "similarity_matrix": similarity_matrix,
    }), 200


# ---------------------------------------------------------------------------
# Entry point (local dev) / Vercel handler
# ---------------------------------------------------------------------------

# Vercel needs 'app' as the WSGI application
handler = app

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
