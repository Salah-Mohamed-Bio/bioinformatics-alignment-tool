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

from alignment_logic import global_alignment, local_alignment, get_symbol

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


@app.route("/align/batch", methods=["POST"])
def align_batch():
    """
    Run batch all-vs-all alignment from a multi-FASTA file upload.
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

    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if j < i:
                row.append(matrix[j][i])
            else:
                _, _, score = global_alignment(seqs[i], seqs[j], match, mismatch, gap)
                row.append(score)
        matrix.append(row)

    return jsonify({
        "ids": ids,
        "matrix": matrix,
    }), 200


# ---------------------------------------------------------------------------
# Entry point (local dev) / Vercel handler
# ---------------------------------------------------------------------------

# Vercel needs 'app' as the WSGI application
handler = app

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)