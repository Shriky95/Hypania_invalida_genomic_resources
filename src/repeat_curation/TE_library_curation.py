#!/usr/bin/env python3
"""
Automated Transposable Element Library Curation Pipeline - Stage 1
Heuristic TE Library Curation

This script executes a multi-step filtering strategy to eliminate spurious models,
host housekeeping gene contamination, and structural classification clashes from 
de novo repeat consensus libraries.

Author: Shrikant Vichare, Bio-Med-Chem Doctoral School, University of Lodz and Lodz Institutes of the Polish Academy of Sciences, Banacha 12/16, 90-237 Lodz, Poland.
License: CC-BY-4.0 / MIT
"""

import pandas as pd
import os
import sys
import argparse

def evaluate_row_triage(row, leak_log_list):
    # 1. Normalization and Cleaning of Core String Attributes
    rm_class = str(row['RepModClass']).strip().lower() if pd.notna(row['RepModClass']) else "unknown"
    dt_class = str(row['DeepTE_Class']).strip().lower() if pd.notna(row['DeepTE_Class']) else "unknown"
    pfam = str(row['Pfam_Domain']).strip().lower() if pd.notna(row['Pfam_Domain']) else "na"
    dfam_model = str(row['Dfam_Model']).strip().lower() if pd.notna(row['Dfam_Model']) else "na"
    target_protein = str(row['Heuristic_Target_Protein']).strip().lower() if pd.notna(row['Heuristic_Target_Protein']) else "na"
    cd_hit = str(row['CD-HIT']).strip().upper() if pd.notna(row['CD-HIT']) else "NR"
    rbh_verdict = str(row['RBH_Verdict']).strip().lower() if 'RBH_Verdict' in row and pd.notna(row['RBH_Verdict']) else "na"
    fam_id = str(row['Family_ID']).strip()
    proposed_nom = str(row['Proposed_Nomenclature']).strip() if 'Proposed_Nomenclature' in row else fam_id
    
    # 2. Strict Numeric Conversion Blocks (Type Safety Gates)
    try:
        seq_len = float(row['Sequence_Length']) if pd.notna(row['Sequence_Length']) else 0.0
    except ValueError:
        seq_len = 0.0
        
    try:
        dfam_bit = float(row['Dfam_BitScore']) if pd.notna(row['Dfam_BitScore']) and str(row['Dfam_BitScore']).upper() != 'NA' else 0.0
    except ValueError:
        dfam_bit = 0.0

    try:
        fwd_qcov = float(row['Fwd_QCovHSP_Max']) if pd.notna(row['Fwd_QCovHSP_Max']) and str(row['Fwd_QCovHSP_Max']).upper() != 'NA' else 0.0
    except ValueError:
        fwd_qcov = 0.0

    # STEP 1: Hard Baseline Filter Pass
    if seq_len < 80.0 or cd_hit == "R":
        return "PURGE_BASELINE"

    # STEP 2: Hard-Gate Manual Curation Discard Verdicts
    if "leak" in rbh_verdict or "purge" in rbh_verdict or "host" in rbh_verdict:
        leak_log_list.append({"id": fam_id, "name": proposed_nom, "domain": f"RBH_Verdict_Flag: {rbh_verdict}"})
        return "PURGE_HOST_CONTAMINATION"

    # Compile Segregated Context Layers to Prevent False Machinery Rescues
    global_context = f"{rm_class} {dt_class} {pfam} {dfam_model} {target_protein}"
    functional_evidence = f"{pfam} {dfam_model} {target_protein}"

    # STEP 3: Expanded Invertebrate Housekeeping & Metabolic Contamination Matrix
    host_keywords = [
        'actin', 'tubulin', 'myosin', 'keratin', 'collagen', 'histone', 'ribosomal', 
        'rrna', 'lsu-rrna', 'ssu-rrna', 'ubiquitin', 'gapdh', 'cytochrome', 'enolase',
        'atp synthase', 'elongation factor', 'translation factor', 'phosphatase', 
        'dehydrogenase', 'nadar', 'pkinase', 'kinase', 'ank_2', 'ankyrin', 'beta-prop', 
        'fbx', 'ift122', 'bloc1', 'clr5', 'crlf3', 'ctd9', 'cyclin', 'ef-hand', 
        'extensin', 'fam81', 'forkhead', 'glyco_transf', 'gst_n', 'hiran', 'homeobox', 
        'inpase', 'irk', 'potassium channel', 'jacalin', 'ldl_recept', 'lrr_8', 'macro', 
        'mob_synth', 'p2rx7', 'purinoceptor', 'polymerase i', 'polymerase ii', 'polymerase iii', 'innexin',
        'peptidase', 'protease', 'subtilisin', 'hydrolase', 'transferase', 'synth synthase', 
        'carboxypeptidase', 'fibrinogen', 'transporter', 'monocarboxylate', 'tropomyosin', 'peroxidase'
    ]
    
    te_rescue_keywords = [
        'transposase', 'reverse transcriptase', 'integrase', 'gag', 'pol', 'env', 
        'thap', 'cenpb', 'cenp-b', 'tesmin', 'mariner', 'hat', 'gypsy', 'copia', 
        'pogo', 'brkdbd', 'dde_', 'tnp_', 'duf1759', 'duf8422', 'exo_endo_phos', 
        'dimer_tnp_hat', 'zbed', 'giy-yig', 'harbi1', 'rve', 'rvt_1', 'rvt_2'
    ]

    has_host_signal = any(hkw in global_context for hkw in host_keywords)
    has_te_machinery = any(ekw in functional_evidence for ekw in te_rescue_keywords)
    
    # Custom Accessory Helicase Exception
    is_helicase_rescue = ('dead' in pfam or 'helicase' in target_protein or 'recq' in target_protein or 'werner' in target_protein) and \
                         ('academ' in global_context or 'harbinger' in global_context or 'pif' in global_context)

    if has_host_signal and not has_te_machinery and not is_helicase_rescue:
        leak_log_list.append({"id": fam_id, "name": proposed_nom, "domain": pfam})
        return "PURGE_HOST_CONTAMINATION"

    # STEP 4: Quality Homology Filters constraining pass conditions to Canonical TE Domains
    canonical_te_domains = [
        'rvt', 'rve', 'integrase', 'transposase', 'gag', 'pol', 'env', 'thap', 
        'cenpb', 'mariner', 'hat', 'gypsy', 'copia', 'pogo', 'brkdbd', 'dde', 'tnp',
        'endonuclease', 'exonuclease', 'exo_endo_phos', 'duf1759', 'duf8422'
    ]
    is_canonical_te_pfam = any(domain in pfam for domain in canonical_te_domains)

    # Core Logic Shift: Non-TE signatures cannot skip checks without >25.0 bitscore
    has_homology_validation = (pfam != "na" and is_canonical_te_pfam) or (dfam_bit >= 25.0)

    broad_te_protein_terms = [
        'reverse transcriptase', 'transposase', 'integrase', 'polyprotein', 
        'piggybac', 'endonuclease', 'exonuclease', 'mobile element', 'transposable element'
    ]
    has_strong_protein_rescue = any(term in target_protein for term in broad_te_protein_terms) and (fwd_qcov >= 50.0)
    has_borderline_dfam_rescue = (dfam_bit >= 20.0)

    if not has_homology_validation and not is_helicase_rescue and not has_strong_protein_rescue and not has_borderline_dfam_rescue:
        if dfam_bit > 0.0 and dfam_bit < 25.0:
            return "PURGE_SPURIOUS_HOMOLOGY"
        
        # De Novo Structural Component Pipeline Concordance Validation Layer
        is_rm_class_dna = any(x in rm_class for x in ['dna', 'hat', 'mariner', 'mu', 'pif', 'harbinger', 'helitron'])
        is_dt_class_dna = any(x in dt_class for x in ['classii', 'dna', 'mite', 'helitron'])
        
        is_rm_class_rna = any(x in rm_class for x in ['ltr', 'gypsy', 'copia', 'pao', 'bel', 'line', 'sine', 'penelope', 'ple'])
        is_dt_class_rna = any(x in dt_class for x in ['classi', 'ltr', 'line', 'sine', 'ple'])
        
        concordant_dna = is_rm_class_dna and is_dt_class_dna
        concordant_rna = is_rm_class_rna and is_dt_class_rna
        
        if not (concordant_dna or concordant_rna or has_te_machinery):
            return "PURGE_STRUCTURAL_CLASH"

    return "RETAIN"

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Empirical Transposable Element Triage and Purging Engine")
    parser.add_argument("-i", "--input", default="Hypania_invalida_master_curation.tsv", help="Path to input master curation TSV file")
    parser.add_argument("-o", "--output", default="retained_elements.tsv", help="Path to output retained elements TSV file")
    parser.add_argument("-l", "--leak_log", default="structural_leaks_purged.txt", help="Path to output text log for purged assets")
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"[-] Operational Error: Input data file layout missing -> '{args.input}'")
        sys.exit(1)

    print(f"[+] Parsing master curation dataset layout: {args.input}")
    df = pd.read_csv(args.input, sep='\t')
    initial_rows = len(df)
    
    print("[+] Evaluating records through the adaptive heuristic triage pipeline matrix...")
    escaped_leaks = []
    triage_decisions = df.apply(lambda r: evaluate_row_triage(r, escaped_leaks), axis=1)
    
    print("\n[=] Triage Step Metric Results Summary:")
    print("-" * 50)
    print(triage_decisions.value_counts().to_string())
    print("-" * 50)
    
    print(f"\n[!] Identified {len(escaped_leaks)} cellular housekeeping or structural leakage models.")
    with open(args.leak_log, "w") as lf:
        for leak in escaped_leaks:
            lf.write(f"{leak['name']}\n")
            
    retained_df = df[triage_decisions == "RETAIN"].copy()
    
    # Strip any older nomenclature metrics to allow fresh recalculation sweeps
    cols_to_drop = [c for c in ['Proposed_Nomenclature', 'Curation_Remark'] if c in retained_df.columns]
    retained_df.drop(columns=cols_to_drop, errors='ignore', inplace=True)
    
    retained_df.to_csv(args.output, sep='\t', index=False)
    print(f"\n[✔] Stage 1 Complete. Clean framework saved to: '{args.output}' ({len(retained_df)} of {initial_rows} rows retained).\n")

if __name__ == "__main__":
    main()