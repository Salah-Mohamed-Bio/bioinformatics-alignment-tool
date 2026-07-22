from flask import Flask, request, jsonify
from flask_cors import CORS
from Bio import SeqIO
import io
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from alignment_logic import global_alignment

app = Flask(__name__)
CORS(app)  # السماح بالاتصال الخارجي

@app.route("/", methods=["GET"])
def home():
    return jsonify({"status": "Flask Bioinformatics API is running successfully!"})

@app.route("/api/batch-align", methods=["POST"])
def batch_align():
    try:
        match = int(request.form.get("match", 5))
        mismatch = int(request.form.get("mismatch", -1))
        gap = int(request.form.get("gap", -3))
        
        if "file" not in request.files:
            return jsonify({"error": "No FASTA file provided"}), 400
            
        file = request.files["file"]
        file_content = file.read().decode("utf-8")
        
        stringio = io.StringIO(file_content)
        sequences = list(SeqIO.parse(stringio, "fasta"))
        
        if len(sequences) < 2:
            return jsonify({"error": "At least 2 sequences are required for batch alignment."}), 400
            
        n = len(sequences)
        matrix = np.zeros((n, n))
        seq_ids = [s.id for s in sequences]
        
        for i in range(n):
            for j in range(i, n):
                _, _, score = global_alignment(str(sequences[i].seq), str(sequences[j].seq), match, mismatch, gap)
                matrix[i, j] = matrix[j, i] = score
                
        return jsonify({
            "success": True,
            "ids": seq_ids,
            "matrix": matrix.tolist()
        })
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True)