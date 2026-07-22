"""
Flask API server for sequence alignment tool.
Wraps alignment_logic.py with REST endpoints for pairwise and batch alignment.
Run locally with: python server.py
"""
import io
import re
import os
import sys

# Ensure alignment_logic.py is importable from project root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from flask import Flask, request, jsonify, render_template
from alignment_logic import global_alignment, local_alignment, get_symbol

app = Flask(__name__)

# Increase max upload size to 16MB
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024

# Maximum sequence length to prevent memory issues on free tier
MAX_SEQ_LENGTH = 5000


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_fasta_text(text):
    """Parse a FASTA string and return the sequence."""
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
    """Run global or local alignment and return results."""
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
    return {
        "align1": align1,
        "align2": align2,
        "symbol": symbol,
        "score": score,
        "mode": mode,
    }


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/", methods=["GET"])
def index():
    """Serve the web frontend."""
    public_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "public", "index.html")
    with open(public_path, "r", encoding="utf-8") as f:
        return f.read()


@app.route("/api/health", methods=["GET"])
def health():
    """Health check endpoint."""
    return jsonify({"status": "ok"})


@app.route("/api/align", methods=["POST"])
def align():
    """
    Run sequence alignment from JSON payload.

    Expected JSON body:
    {
        "seq1": "ACGT...",
        "seq2": "TGCA...",
        "match": 5,
        "mismatch": -1,
        "gap": -3,
        "mode": "global" | "local"
    }
    """
    data = request.get_json(silent=True)
    if not data:
        return jsonify({"error": "Request body must be valid JSON."}), 400

    seq1 = data.get("seq1", "").strip().upper()
    seq2 = data.get("seq2", "").strip().upper()
    match = data.get("match", 5)
    mismatch = data.get("mismatch", -1)
    gap = data.get("gap", -3)
    mode = data.get("mode", "global")

    try:
        result = run_alignment(seq1, seq2, match, mismatch, gap, mode)
        return jsonify(result), 200
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except MemoryError:
        return jsonify({"error": "Not enough memory for this alignment. Please use shorter sequences."}), 400


@app.route("/api/align/fasta", methods=["POST"])
def align_fasta():
    """
    Run sequence alignment from uploaded FASTA files.

    Multipart form fields:
        file1 - first sequence FASTA file
        file2 - second sequence FASTA file
        match - match score (default: 5)
        mismatch - mismatch score (default: -1)
        gap - gap penalty (default: -3)
        mode - "global" or "local" (default: "global")
    """
    if "file1" not in request.files or "file2" not in request.files:
        return jsonify({"error": "Both 'file1' and 'file2' are required."}), 400

    try:
        seq1 = parse_fasta_file(request.files["file1"])
        seq2 = parse_fasta_file(request.files["file2"])
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    match = request.form.get("match", 5, type=int)
    mismatch = request.form.get("mismatch", -1, type=int)
    gap = request.form.get("gap", -3, type=int)
    mode = request.form.get("mode", "global")

    try:
        result = run_alignment(seq1, seq2, match, mismatch, gap, mode)
        return jsonify(result), 200
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except MemoryError:
        return jsonify({"error": "Not enough memory for this alignment. Please use shorter sequences."}), 400


@app.route("/api/align/batch", methods=["POST"])
def align_batch():
    """
    Run batch all-vs-all alignment from a multi-FASTA file upload.

    Multipart form fields:
        file - multi-FASTA file with multiple sequences (2-20)
        match - match score (default: 5)
        mismatch - mismatch score (default: -1)
        gap - gap penalty (default: -3)

    Returns JSON:
    {
        "ids": ["seq1", "seq2", ...],
        "matrix": [[score, ...], ...]
    }
    """
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

    match = request.form.get("match", 5, type=int)
    mismatch = request.form.get("mismatch", -1, type=int)
    gap = request.form.get("gap", -3, type=int)

    ids = [r[0] for r in records]
    seqs = [r[1] for r in records]
    n = len(seqs)

    # Build similarity matrix (all-vs-all)
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if j < i:
                row.append(matrix[j][i])  # mirror symmetric
            else:
                _, _, score = global_alignment(seqs[i], seqs[j], match, mismatch, gap)
                row.append(score)
        matrix.append(row)

    return jsonify({
        "ids": ids,
        "matrix": matrix,
    }), 200


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
