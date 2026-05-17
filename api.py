from fastapi import FastAPI, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
import pandas as pd
import tempfile
import os
from vaccine_builder import (
    load_epitopes,
    score_candidates,
    select_targets,
    assemble_polyepitope,
    build_mrna_construct
)

app = FastAPI()

# This allows the Lovable frontend to talk to this backend
# without getting blocked by browser security rules
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def root():
    return {"status": "OpenVAXX API is running"}

@app.post("/design")
async def design_vaccine(file: UploadFile = File(...)):

    # Save the uploaded file temporarily so we can read it
    with tempfile.NamedTemporaryFile(delete=False, suffix=".tsv") as tmp:
        contents = await file.read()
        tmp.write(contents)
        tmp_path = tmp.name

    try:
        # Run the full pipeline on the uploaded file
        df = load_epitopes(tmp_path)
        df_scored = score_candidates(df)
        targets = select_targets(df_scored, n=10)
        polyepitope = assemble_polyepitope(targets)
        construct = build_mrna_construct(polyepitope)

        # Build the response with everything the frontend needs
        epitope_list = []
        for _, row in targets.iterrows():
            epitope_list.append({
                "gene": str(row.get("Gene", "")),
                "mutation": str(row.get("AA Change", "")),
                "peptide": str(row.get("Best Peptide", "")),
                "ic50": round(float(row.get("IC50 MT", 0)), 2),
                "fold_change": round(float(row.get("fold_change", 0)), 2),
                "dna_vaf": round(float(row.get("DNA VAF", 0) or 0), 3),
                "priority_score": round(float(row.get("priority_score", 0)), 3),
            })

        return {
            "status": "success",
            "num_candidates_analyzed": len(df),
            "num_targets_selected": len(targets),
            "polyepitope": polyepitope,
            "construct_length": len(construct),
            "construct": construct,
            "epitopes": epitope_list,
        }

    except Exception as e:
        return {"status": "error", "message": str(e)}

    finally:
        # Clean up the temporary file
        os.unlink(tmp_path)
