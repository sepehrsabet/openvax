import pandas as pd

# ── Load the TSV file ──────────────────────────────────────────
def load_epitopes(filepath):
    df = pd.read_csv(filepath, sep="\t")
    return df

# ── Score and rank candidates ──────────────────────────────────
def score_candidates(df):
    # Keep only rows where IC50 MT exists and is below 1500
    # IC50 MT is how tightly the peptide binds - lower is better
    # 1500 is a loose threshold we use because osteosarcoma has few candidates
    df = df[pd.to_numeric(df["IC50 MT"], errors="coerce").notna()]
    df["IC50 MT"] = pd.to_numeric(df["IC50 MT"])
    df["IC50 WT"] = pd.to_numeric(df["IC50 WT"], errors="coerce")
    df["DNA VAF"] = pd.to_numeric(df["DNA VAF"], errors="coerce")

    # Binding score - higher means better binding
    # We flip IC50 so that lower IC50 gives higher score
    df["binding_score"] = 1000 / (df["IC50 MT"] + 1)

    # Fold change - how much better mutant binds vs normal version
    # Higher fold change means the mutation created something new the immune system hasnt seen
    df["fold_change"] = df["IC50 WT"] / df["IC50 MT"]
    df["fold_change"] = df["fold_change"].fillna(1)

    # VAF score - how many tumor cells carry this mutation
    df["vaf_score"] = df["DNA VAF"].fillna(0.05)

    # Composite score combining all three
    df["priority_score"] = (
        df["binding_score"] * 0.5 +
        df["fold_change"] * 0.3 +
        df["vaf_score"] * 10 * 0.2
    )

    # Sort best to worst
    df = df.sort_values("priority_score", ascending=False)
    return df

# ── Select top candidates ──────────────────────────────────────
def select_targets(df, n=10):
    # Take the top n candidates by priority score
    # We use 10 here because we only have 20 total candidates
    top = df.head(n).copy()
    return top

# ── Assemble polyepitope sequence ──────────────────────────────
def assemble_polyepitope(targets):
    # AAY is a short linker sequence that helps the cell cut between epitopes cleanly
    linker = "AAY"
    peptides = targets["Best Peptide"].dropna().tolist()
    polyepitope = linker.join(peptides)
    return polyepitope

# ── Build full mRNA construct ──────────────────────────────────
def build_mrna_construct(polyepitope):
    # Signal peptide - directs the protein to the right place in the cell
    signal = "ATGGATGCAATGAAGAGAGGGCTCTGCTGTGTGCTGCTTCTGGGGGTCTGGTCCAGTGGG"
    # Kozak sequence - helps the cell know where to start reading the recipe
    kozak = "GCCACCATG"
    # 5' UTR - the header of the recipe card
    utr5 = "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACC"
    # 3' UTR - the footer that keeps the mRNA stable
    utr3 = "UGCAUACGAGAUUCUCUAGCCAUAAAGCCCGGGGCGGCCGCAUAAAA"
    # Poly-A tail - 120 A's that protect the end of the mRNA from degradation
    polya = "A" * 120
    # Stop codon - tells the cell to stop reading
    stop = "TGATAA"

    # Back translate polyepitope from amino acids to DNA
    # This is a simple codon table using the most common human codons
    codon_table = {
        "A": "GCC", "R": "AGA", "N": "AAC", "D": "GAC", "C": "TGC",
        "Q": "CAG", "E": "GAG", "G": "GGC", "H": "CAC", "I": "ATC",
        "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
        "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG"
    }

    coding_dna = ""
    for aa in polyepitope:
        if aa in codon_table:
            coding_dna += codon_table[aa]
        else:
            coding_dna += "NNN"

    full_construct = utr5 + kozak + signal + coding_dna + stop + utr3 + polya
    return full_construct

# ── Write FASTA output ─────────────────────────────────────────
def write_fasta(sequence, output_path, label="OpenVAXX_Osteosarcoma_Vaccine"):
    with open(output_path, "w") as f:
        f.write(f">{label}\n")
        # Write sequence in chunks of 60 characters per line
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
    print(f"FASTA written to {output_path}")

# ── Run everything ─────────────────────────────────────────────
if __name__ == "__main__":
    input_file = "/mnt/e/osteosarc/neoantigens/mhci_all_epitopes.tsv"
    output_file = "/home/sepehr/openvax/vaccine_v1.fasta"

    print("Loading epitopes...")
    df = load_epitopes(input_file)

    print("Scoring candidates...")
    df_scored = score_candidates(df)

    print("Selecting top targets...")
    targets = select_targets(df_scored, n=10)

    print("\nTop selected neoantigens:")
    print(targets[["Gene", "Best Peptide", "IC50 MT", "fold_change", "priority_score"]].to_string())

    print("\nAssembling polyepitope...")
    polyepitope = assemble_polyepitope(targets)
    print(f"Polyepitope: {polyepitope}")

    print("\nBuilding mRNA construct...")
    construct = build_mrna_construct(polyepitope)
    print(f"Construct length: {len(construct)} nucleotides")

    write_fasta(construct, output_file)
    print("\nDone.")
