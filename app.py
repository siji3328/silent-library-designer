# -*- coding: utf-8 -*-
"""Silent Library Designer - synonymous DNA library oligomer designer.

A single-tool web app: paste a WT coding sequence and its amino-acid sequence, get back a
degenerate (IUPAC) synonymous library split into PCA-ready oligomers. The design logic
lives in tools/; this module is only routing and form handling.
"""

import os

from flask import Flask, Response, render_template, request

from tools import silent_library
from tools.common import pca_rows_to_csv, pca_rows_to_fasta

app = Flask(__name__)

VERSION = "1.0"


@app.route("/", methods=["GET", "POST"])
def designer():
    form = {
        "dna_sequence": silent_library.DEFAULT_DNA_SEQUENCE,
        "amino_acid_sequence": silent_library.DEFAULT_AA_SEQUENCE,
        "oligo_length": 89,
        "overlap": 20,
        "auto_split": True,
    }
    result = None
    error = None
    warning = None

    if request.method == "POST":
        form["dna_sequence"] = request.form.get("dna_sequence", "").strip().upper()
        form["amino_acid_sequence"] = request.form.get("amino_acid_sequence", "").strip().upper()
        form["oligo_length"] = request.form.get("oligo_length", "89")
        form["overlap"] = request.form.get("overlap", "20")
        form["auto_split"] = request.form.get("auto_split") == "on"

        dna_sequence = "".join(form["dna_sequence"].split())
        amino_acid_sequence = "".join(form["amino_acid_sequence"].split())

        auto_split = form["auto_split"]
        try:
            oligo_length = int(form["oligo_length"])
            overlap = int(form["overlap"])
        except ValueError:
            oligo_length = overlap = None
            error = "올리고 길이와 오버랩 길이는 정수로 입력해주세요."

        candidate_key = request.form.get("candidate") or None

        if not error:
            result, error, warning = silent_library.run(
                dna_sequence, amino_acid_sequence, oligo_length, overlap, auto_split, candidate_key
            )

        action = request.form.get("action", "run")
        if result and action in ("download_csv", "download_fasta"):
            if action == "download_csv":
                body = pca_rows_to_csv(result["oligos"])
                return Response(
                    body, mimetype="text/csv",
                    headers={"Content-Disposition": "attachment; filename=silent_library_oligos.csv"},
                )
            body = pca_rows_to_fasta(result["oligos"])
            return Response(
                body, mimetype="text/x-fasta",
                headers={"Content-Disposition": "attachment; filename=silent_library_oligos.fasta"},
            )

    return render_template(
        "index.html", form=form, result=result, error=error, warning=warning, version=VERSION
    )


if __name__ == "__main__":
    # Port 7860 matches the Docker/Hugging Face Spaces setup, so local and deployed runs agree.
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 7860)), debug=True)
