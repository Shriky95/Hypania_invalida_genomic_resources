#!/usr/bin/env python3
"""
Automated Transposable Element Library Curation Pipeline - Stage 2
Empirical Wicker Nomenclature

This script parses curated transposable element families and maps them into 
standardized Wicker taxonomy codes and structural suffixes based on model contiguity.

Author: Shrikant Vichare, Bio-Med-Chem Doctoral School, University of Lodz and Lodz Institutes of the Polish Academy of Sciences, Banacha 12/16, 90-237 Lodz, Poland.
License: CC-BY-4.0 / MIT
"""

import pandas as pd
import re
import os
import sys
import argparse

def assign_wicker_profile(row, species_prefix):
    fam_id = str(row['Family_ID']).strip()
    rm_class = str(row['RepModClass']).strip() if pd.notna(row['RepModClass']) else "Unknown"
    dt_class = str(row['DeepTE_Class']).strip() if pd.notna(row['DeepTE_Class']) else "unknown"
    pfam = str(row['Pfam_Domain']).strip() if pd.notna(row['Pfam_Domain']) else "NA"
    dfam_model = str(row['Dfam_Model']).strip() if pd.notna(row['Dfam_Model']) else "NA"
    target_protein = str(row['Heuristic_Target_Protein']).strip() if pd.notna(row['Heuristic_Target_Protein']) else "NA"
    
    try:
        seq_len = float(row['Sequence_Length']) if pd.notna(row['Sequence_Length']) else 0.0
    except ValueError:
        seq_len = 0.0
        
    try:
        dfam_bit = float(row['Dfam_BitScore']) if pd.notna(row['Dfam_BitScore']) and str(row['Dfam_BitScore']).upper() != 'NA' else 0.0
    except ValueError:
        dfam_bit = 0.0

    id_match = re.search(r'rnd-(\d+)_family-(\d+)', fam_id)
    suffix_index = f"_{id_match.group(1)}_{id_match.group(2)}" if id_match else "_0_0"
    
    rm_low, dt_low, pf_low, df_low, tp_low = rm_class.lower(), dt_class.lower(), pfam.lower(), dfam_model.lower(), target_protein.lower()
    unified_context = f"{rm_low} {pf_low} {df_low} {tp_low} {dt_low}"

    wicker_code, superfamily = "UXX", "Unknown"
    evidence_trail = "Classification resolved via global sequence component profiling."

    # -------------------------------------------------------------------------
    #               HIERARCHICAL TAXONOMY RESOLUTION MATRIX
    # -------------------------------------------------------------------------
    
    # Tier 1: Prioritize explicit, high-confidence RepeatMasker definitions
    if "mariner" in rm_low or "tc1" in rm_low or "tc2" in rm_low or "tcmar" in rm_low or "pogo" in rm_low:
        wicker_code, superfamily = "DTT", "TcMariner"
        evidence_trail = "Class II DNA Tc1-Mariner lineage enforced via explicit RepeatMasker taxonomy parameters."
    elif "hat" in rm_low:
        wicker_code, superfamily = "DTA", "hAT"
        evidence_trail = "Class II DNA hAT superfamily classification verified via explicit RepeatMasker parameters."
    elif "mutator" in rm_low or "mule" in rm_low or "mu-" in rm_low:
        wicker_code, superfamily = "DTM", "Mutator"
        evidence_trail = "Class II DNA Mutator element verified via explicit RepeatMasker validation markers."
    elif "cacta" in rm_low or "enspm" in rm_low or "en-spm" in rm_low or "chapaev" in rm_low or "cmc" in rm_low:
        wicker_code, superfamily = "DTC", "CMC_Chapaev" if ("chapaev" in rm_low or "chapaev" in pf_low) else "CACTA"
        evidence_trail = "Class II DNA CACTA-Chapaev array verified via explicit RepeatMasker taxonomy parameters."
    elif "harbinger" in rm_low or "pif-" in rm_low:
        wicker_code, superfamily = "DTH", "Harbinger"
        evidence_trail = "Class II DNA PIF-Harbinger structural profile verified via primary RepeatMasker markers."
    elif "academ" in rm_low:
        wicker_code, superfamily = "DXX", "Academ"
        evidence_trail = "Class II DNA Academ lineage enforced via primary RepeatMasker markers."
    elif "helitron" in rm_low or "rc/" in rm_low:
        wicker_code, superfamily = "DHH", "Helitron"
        evidence_trail = "Class III Helitron rolling-circle framework validated via RepeatMasker profile flags."
    elif "maverick" in rm_low or "polinton" in rm_low:
        wicker_code, superfamily = "DMX", "Maverick"
        evidence_trail = "Class II Maverick/Polinton structural block validated via primary RepeatMasker definitions."
    elif "gypsy" in rm_low:
        wicker_code, superfamily = "RLG", "Gypsy"
        evidence_trail = "Class I LTR retrotransposon Gypsy framework confirmed via primary RepeatMasker definitions."
    elif "copia" in rm_low:
        wicker_code, superfamily = "RLR", "Copia"
        evidence_trail = "Class I LTR retrotransposon Copia framework confirmed via primary RepeatMasker definitions."
    elif "pao" in rm_low or "bel" in rm_low:
        wicker_code, superfamily = "RLB", "BEL"
        evidence_trail = "Class I LTR retrotransposon BEL/Pao framework confirmed via primary RepeatMasker definitions."
    elif "cr1" in rm_low:
        wicker_code, superfamily = "RIL", "CR1"
        evidence_trail = "Class I Non-LTR LINE CR1 framework validated via primary RepeatMasker database profiles."
    elif "l2" in rm_low:
        wicker_code, superfamily = "RIL", "L2"
        evidence_trail = "Class I Non-LTR LINE L2 framework validated via primary RepeatMasker database profiles."
    elif "l1" in rm_low:
        wicker_code, superfamily = "RIL", "L1"
        evidence_trail = "Class I Non-LTR LINE L1 framework validated via primary RepeatMasker database profiles."
    elif "rte" in rm_low or "bovb" in rm_low:
        wicker_code, superfamily = "RIR", "RTE"
        evidence_trail = "Class I Non-LTR LINE RTE-BovB framework validated via primary RepeatMasker parameters."
    elif "jockey" in rm_low:
        wicker_code, superfamily = "RIX", "Jockey"
        evidence_trail = "Class I Non-LTR LINE Jockey framework validated via primary RepeatMasker parameters."
    elif "sine" in rm_low:
        wicker_code, superfamily = "RSX", "SINE"
        evidence_trail = "Class I SINE non-autonomous retrotransposon tracking confirmed via primary parameters."
    elif "penelope" in rm_low or "ple" in rm_low:
        wicker_code, superfamily = "RIX", "Penelope"
        evidence_trail = "Penelope-like element (PLE) tracking confirmed via primary structural parameters."
    elif "line" in rm_low:
        wicker_code, superfamily = "RIX", "LINE"
        evidence_trail = "Class I Non-LTR LINE general distribution validated via primary RepeatMasker profiles."
    elif "ltr" in rm_low:
        wicker_code, superfamily = "RLX", "LTR"
        evidence_trail = "Class I LTR general distribution validated via primary RepeatMasker profiles."

    # Tier 2: Fallback to high-confidence Functional Domains if classifications are broad
    if wicker_code == "UXX" or superfamily == "Unknown":
        if "polinton" in unified_context or "maverick" in unified_context: wicker_code, superfamily = "DMX", "Maverick"
        elif "academ" in unified_context: wicker_code, superfamily = "DXX", "Academ"
        elif "harbinger" in unified_context or "pif-" in unified_context or "harbi1" in unified_context or "dde_tnp_4" in unified_context or "isl2eu" in unified_context: wicker_code, superfamily = "DTH", "Harbinger"
        elif "cacta" in unified_context or "chapaev" in unified_context or "cmc" in unified_context or "enspm" in unified_context or "scaaqu" in unified_context: wicker_code, superfamily = "DTC", "CMC_Chapaev" if ("chapaev" in unified_context or "scaaqu" in unified_context) else "CACTA"
        elif "sine" in unified_context and "sinensis" not in unified_context: wicker_code, superfamily = "RSX", "SINE"
        elif "penelope" in unified_context or "ple" in unified_context or "duf8422" in unified_context: wicker_code, superfamily = "RIX", "Penelope"
        elif "gypsy" in unified_context or "rvt_1" in unified_context: wicker_code, superfamily = "RLG", "Gypsy"
        elif "copia" in unified_context or "rvt_2" in unified_context: wicker_code, superfamily = "RLR", "Copia"
        elif any(x in unified_context for x in ['pao', 'bel']) or "duf1759" in unified_context: wicker_code, superfamily = "RLB", "BEL"
        elif "l1" in unified_context: wicker_code, superfamily = "RIL", "L1"
        elif "l2" in unified_context or "exo_endo_phos" in unified_context: wicker_code, superfamily = "RIL", "L2"
        elif "cr1" in unified_context: wicker_code, superfamily = "RIL", "CR1"
        elif "rte" in unified_context: wicker_code, superfamily = "RIR", "RTE"
        elif "jockey" in unified_context: wicker_code, superfamily = "RIX", "Jockey"
        elif "line" in unified_context: wicker_code, superfamily = "RIX", "LINE"
        elif any(x in unified_context for x in ['mariner', 'tc1', 'tc2', 'tc3', 'tcmar', 'pogo']) or any(y in unified_context for y in ['brkdbd', 'cenp-b', 'hth_tnp_tc3_2', 'dde_1', 'dde_3']): wicker_code, superfamily = "DTT", "TcMariner"
        elif "hat" in unified_context or any(z in unified_context for z in ['dimer_tnp_hat', 'zbed', 'zf-bed']): wicker_code, superfamily = "DTA", "hAT"
        elif "mutator" in unified_context or "mule" in unified_context: wicker_code, superfamily = "DTM", "Mutator"
        elif "helitron" in unified_context: wicker_code, superfamily = "DHH", "Helitron"
        elif "ltr" in dt_low: wicker_code, superfamily = "RLX", "LTR"
        elif "classii" in dt_low or "dna" in dt_low: wicker_code, superfamily = "DXX", "DNA"
        elif "classi" in dt_low: wicker_code, superfamily = "RXX", "Retrotransposon"

    # Tier 3: Broad Base Catch Fallbacks if still unmapped
    if wicker_code == "UXX":
        if "dna" in rm_low: wicker_code, superfamily = "DXX", "DNA"
        elif "ltr" in rm_low: wicker_code, superfamily = "RLX", "LTR"
        elif "line" in rm_low: wicker_code, superfamily = "RIX", "LINE"
        elif "sine" in rm_low: wicker_code, superfamily = "RSX", "SINE"

    # Structural Nomenclature Suffix Formatting Control Loop
    has_coding_domain = any(dom in pf_low for dom in ['dde', 'tnp', 'rvt', 'rve', 'integrase', 'exo_endo_phos', 'brkdbd', 'dimer_tnp_hat', 'dead'])
    if seq_len > 2000.0 or has_coding_domain or wicker_code in ["UXX", "DMX"]:
        struct_suffix = ""
    elif wicker_code.startswith("D"):
        struct_suffix = "_MITE" if seq_len <= 600.0 else "_nMITE"
    else:
        struct_suffix = "_Fragment"

    proposed_nomenclature = f"{species_prefix}_{wicker_code}_{superfamily}{struct_suffix}{suffix_index}".replace("__", "_")
    curation_remark = f"{evidence_trail} Length: {int(seq_len)} bp. Reference Dfam BitScore: {dfam_bit}."
    return proposed_nomenclature, curation_remark

def main():
    parser = argparse.ArgumentParser(description="Stage 2: Empirical Wicker Structural Nomenclature Engine")
    parser.add_argument("-i", "--input", default="retained_elements.tsv", help="Path to whitelisted input TSV from Stage 1")
    parser.add_argument("-o", "--output", default="Hypania_invalida_final_curation.tsv", help="Path to output taxonomy finalized TSV")
    parser.add_argument("-p", "--prefix", default="Hinv", help="Organism specific string identifier prefix (e.g., Hinv)")
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"[-] Operational Error: Input file missing -> '{args.input}'. Ensure Stage 1 has executed successfully.")
        sys.exit(1)

    df = pd.read_csv(args.input, sep='\t')
    print(f"[+] Stage 2: Deploy Wicker taxonomy matrix over {len(df)} elements...")
    
    taxonomy_results = df.apply(lambda row: assign_wicker_profile(row, args.prefix), axis=1)
    df['Proposed_Nomenclature'] = [r[0] for r in taxonomy_results]
    df['Curation_Remark'] = [r[1] for r in taxonomy_results]
    
    print("\n[=] Final Wicker Class Profiles:")
    print("-" * 50)
    print(df['Proposed_Nomenclature'].str.extract(r'_([A-Z]{3})_')[0].value_counts().to_string())
    print("-" * 50)
    
    df.to_csv(args.output, sep='\t', index=False)
    print(f"[✔] Stage 2 Complete. Curated table saved to: '{args.output}'\n")

if __name__ == "__main__":
    main()