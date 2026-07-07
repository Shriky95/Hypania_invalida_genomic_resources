#!/usr/bin/env python3
"""
Automated Transposable Element Library Curation Pipeline - Pre-Curation Framework
Compiles Steps 1-9

This pipeline transitions from an uncurated pool of raw ab initio consensus models 
down to a structurally resolved, mathematically validated 16-column master curation matrix.

Author: Shrikant Vichare, Bio-Med-Chem Doctoral School, University of Lodz and 
        Lodz Institutes of the Polish Academy of Sciences, Banacha 12/16, 90-237 Lodz, Poland.
License: CC-BY-4.0 / MIT
"""

import os
import sys
import re
import argparse
import pandas as pd
import numpy as np

# Canonical Genetic Code Mapping Dictionary
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def reverse_complement(seq):
    pairs = str.maketrans("ATGCNatgcn", "TACGNtacgn")
    return seq.translate(pairs)[::-1]

def translate_frame(seq, frame_offset):
    peptides = []
    for i in range(frame_offset, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        peptides.append(GENETIC_CODE.get(codon, 'X'))
    return "".join(peptides)

def extract_orfs(protein_seq, min_len):
    raw_orfs = protein_seq.split('*')
    return [orf for orf in raw_orfs if len(orf) >= min_len and 'X' not in orf]

def parse_fasta(fasta_path):
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line.split()[0].replace(">", "")
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences

# =============================================================================
#                       PIPELINE STAGE EXECUTIONS
# =============================================================================

def execute_stage1_translation(input_fasta, output_faa, min_orf_len):
    print(f"[+] Launching Stage 1: Six-Frame In Silico Translatome Profiler on '{input_fasta}'")
    raw_library = parse_fasta(input_fasta)
    print(f"    [>] Loaded {len(raw_library)} nucleotide consensus models.")
    
    viable_peptides_count = 0
    with open(output_faa, 'w') as faa_out:
        for te_id, nt_seq in raw_library.items():
            rc_seq = reverse_complement(nt_seq)
            frames = {
                'F1': translate_frame(nt_seq, 0), 'F2': translate_frame(nt_seq, 1), 'F3': translate_frame(nt_seq, 2),
                'R1': translate_frame(rc_seq, 0), 'R2': translate_frame(rc_seq, 1), 'R3': translate_frame(rc_seq, 2)
            }
            for frame_name, aa_seq in frames.items():
                orfs = extract_orfs(aa_seq, min_orf_len)
                for index, orf in enumerate(orfs, 1):
                    header = f">{te_id}_frame_{frame_name}_orf_{index} [length={len(orf)}]\n"
                    faa_out.write(header + orf + "\n")
                    viable_peptides_count += 1
                    
    print(f"    [✔] Stage 1 Complete. Exported {viable_peptides_count} viable peptides to '{output_faa}'.")

def execute_stage3_heuristic_scoring(fwd_hits_raw, matrix_out, targets_txt_out):
    print(f"[+] Launching Stage 3: Chunked Aggregation & Heuristic Scoring (4:3:2:1 Model) on '{fwd_hits_raw}'")
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'qcovhsp', 'stitle']
    
    chunk_iter = pd.read_csv(fwd_hits_raw, sep='\t', header=None, names=cols, on_bad_lines='skip', chunksize=500000)
    agg_list = []
    
    for chunk in chunk_iter:
        chunk_agg = chunk.groupby(['qseqid', 'sseqid']).agg({
            'pident': 'max', 'evalue': 'min', 'bitscore': ['max', 'sum'], 'length': 'sum', 'qcovhsp': 'max', 'stitle': 'first'
        })
        chunk_agg.columns = ['_'.join(col).strip() for col in chunk_agg.columns.values]
        agg_list.append(chunk_agg.reset_index())
        
    df = pd.concat(agg_list).groupby(['qseqid', 'sseqid']).agg({
        'pident_max': 'max', 'evalue_min': 'min', 'bitscore_max': 'max', 'bitscore_sum': 'sum', 
        'length_sum': 'sum', 'qcovhsp_max': 'max', 'stitle_first': 'first'
    }).reset_index()
    
    df['neg_log_evalue'] = -np.log10(df['evalue_min'] + 1e-300)
    df['norm_bit'] = df.groupby('qseqid')['bitscore_sum'].transform(lambda x: x / x.max() if x.max() > 0 else 0)
    df['norm_eval'] = df.groupby('qseqid')['neg_log_evalue'].transform(lambda x: x / x.max() if x.max() > 0 else 0)
    df['norm_cov'] = df['qcovhsp_max'] / 100.0
    df['norm_idt'] = df['pident_max'] / 100.0
    
    # 4:3:2:1 Matrix Formula Assignment
    df['heuristic_score'] = (0.40 * df['norm_bit']) + (0.30 * df['norm_cov']) + (0.20 * df['norm_idt']) + (0.10 * df['norm_eval'])
    
    df['max_h'] = df.groupby('qseqid')['heuristic_score'].transform('max')
    top10_dynamic = df[df['heuristic_score'] >= (df['max_h'] * 0.70)].sort_values(
        by=['qseqid', 'heuristic_score'], ascending=[True, False]
    ).groupby('qseqid').head(10)
    
    top10_dynamic.to_csv(matrix_out, sep='\t', index=False)
    
    all_target_accessions = set()
    for _, row in top10_dynamic.iterrows():
        s = row['sseqid']
        acc = s.split('|')[-2] if '|' in s else s
        all_target_accessions.add(acc)
        
    with open(targets_txt_out, "w") as f:
        for acc in sorted(all_target_accessions):
            f.write(f"{acc}\n")
            
    print(f"    [✔] Stage 3 Complete. Top matrices saved. Extracted {len(all_target_accessions)} unique reference targets.")

def execute_stage4_stream_database(master_ref_fasta, targets_txt, output_reciprocal_faa):
    print(f"[+] Launching Stage 4: Hashed Database Streaming Extraction on '{master_ref_fasta}'")
    with open(targets_txt) as f:
        targets = set(line.strip() for line in f)
        
    with open(master_ref_fasta, "r") as infile, open(output_reciprocal_faa, "w") as outfile:
        write_seq = False
        for line in infile:
            if line.startswith(">"):
                line_id = line.split()[0][1:]
                write_seq = line_id in targets
                if write_seq: outfile.write(line)
            elif write_seq:
                outfile.write(line)
                
    print(f"    [✔] Stage 4 Complete. Reciprocal queries database constructed at '{output_reciprocal_faa}'.")

def execute_stage5_6_rbh_and_clean(fwd_matrix_tsv, rev_hits_raw, collapsed_report_out):
    print("[+] Launching Stages 5 & 6: Bidirectional RBH Parity Integration & Host Contamination Clean-up")
    fwd_cols = ['qseqid', 'sseqid', 'pident_max', 'evalue_min', 'bitscore_max', 'bitscore_sum', 'length_sum', 'qcovhsp_max', 'stitle_first', 'heuristic_score']
    fwd_df = pd.read_csv(fwd_matrix_tsv, sep='\t')[fwd_cols]
    fwd_best = fwd_df.sort_values(by=['qseqid', 'heuristic_score'], ascending=[True, False]).drop_duplicates(subset=['qseqid']).copy()
    fwd_best['clean_sseqid'] = fwd_best['sseqid'].apply(lambda x: x.split('|')[-2] if '|' in x else x)
    
    rev_cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'qcovhsp', 'stitle']
    rev_df = pd.read_csv(rev_hits_raw, sep='\t', header=None, names=rev_cols)
    rev_df['clean_qseqid'] = rev_df['qseqid'].apply(lambda x: x.split('|')[-2] if '|' in x else x)
    
    merged = pd.merge(fwd_best, rev_df, left_on='clean_sseqid', right_on='clean_qseqid', suffixes=('_fwd', '_rev'))
    merged['is_rbh'] = merged['qseqid_fwd'] == merged['sseqid_rev']
    
    # Collapse rows down to unified Family ID resolution
    merged['clean_family_id'] = merged['qseqid_fwd'].astype(str).str.split('_frame_').str[0].str.split('#').str[0].str.strip()
    
    host_keywords = ['actin', 'tubulin', 'ribosomal', 'histone', 'ubiquitin', 'myosin', 'kinase', 'dehydrogenase', 
                     'polymerase', 'synthase', 'helicase', 'chaperone', 'heat shock', 'cytochrome', 'transcription factor', 'ferritin']
    te_exceptions = ['transposase', 'retrotransposon', 'reverse transcriptase', 'integrase', 'pol protein', 'gag protein']
    
    def compute_verdict(row):
        if row['is_rbh']:
            title = str(row['stitle_first']).lower()
            if any(k in title for k in host_keywords) and not any(e in title for e in te_exceptions):
                return "Host Gene Leak (Purge)"
            return "True RBH (Autonomous TE)"
        return "False RBH (Fragment/Paralog)"
        
    merged['RBH_Verdict'] = merged.apply(compute_verdict, axis=1)
    df_sorted = merged.sort_values(by=['clean_family_id', 'heuristic_score'], ascending=[True, False])
    df_unique = df_sorted.drop_duplicates(subset=['clean_family_id'])
    
    cleaned_report = df_unique[['clean_family_id', 'stitle_first', 'bitscore_sum', 'qcovhsp_max', 'heuristic_score', 'RBH_Verdict']].rename(
        columns={'stitle_first': 'Heuristic_Target_Protein', 'bitscore_sum': 'Fwd_BitScore_Sum', 'qcovhsp_max': 'Fwd_QCovHSP_Max', 'heuristic_score': 'Heuristic_Score'}
    )
    cleaned_report.to_csv(collapsed_report_out, sep='\t', index=False)
    print(f"    [✔] Stages 5 & 6 Complete. Condensed coding reports saved to '{collapsed_report_out}'.")

def execute_stage7_8_9_consolidation(master_template_tsv, clean_report_tsv, clstr_file, library_fasta, final_curation_tsv):
    print("[+] Launching Stages 7, 8 & 9: Multilayered Consolidation, Redundancy Mapping & Base Length Verification")
    master_df = pd.read_csv(master_template_tsv, sep='\t').replace('NA', np.nan)
    clean_rbh = pd.read_csv(clean_report_tsv, sep='\t')
    
    master_df['Family_ID'] = master_df['Family_ID'].astype(str).str.strip()
    clean_rbh['clean_family_id'] = clean_rbh['clean_family_id'].astype(str).str.strip()
    
    cols_to_keep = ['Family_ID', 'RepModClass', 'Pfam_Domain', 'Pfam_Evalue', 'DeepTE_Class', 'tRNAscan', 'Dfam_Model', 'Dfam_Evalue', 'Dfam_BitScore']
    master_trimmed = master_df[cols_to_keep]
    
    merged = pd.merge(master_trimmed, clean_rbh, left_on='Family_ID', right_on='clean_family_id', how='left')
    
    def fill_missing_verdicts(row):
        if pd.notna(row['RBH_Verdict']) and str(row['RBH_Verdict']) != 'nan':
            return row['RBH_Verdict']
        t_str, r_str = str(row['tRNAscan']).lower().strip(), str(row['RepModClass']).lower().strip()
        if (t_str != 'none' and t_str != 'nan' and pd.notna(row['tRNAscan'])) or 'trna' in r_str:
            return "Structural RNA (Purge/SINE-link)"
        return "Non-Coding / Novel Candidate"
        
    merged['RBH_Verdict'] = merged.apply(fill_missing_verdicts, axis=1).drop(columns=['clean_family_id'], errors='ignore')
    
    # CD-HIT Cluster Parsing Logic Block
    status_map = {}
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>'): continue
            match = re.search(r'>([^ .]+)', line)
            if match:
                f_id = match.group(1).split('#')[0].strip()
                is_rep = '*' in line
                if f_id not in status_map or is_rep:
                    status_map[f_id] = 'NR' if is_rep else 'R'
                    
    merged['CD-HIT'] = merged['Family_ID'].map(status_map).fillna('NR')
    
    # FASTA Physical Base Pair Metric Assembly Block
    fasta_lengths = {}
    current_family_id, current_seq_accumulator = None, []
    with open(library_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_family_id and current_seq_accumulator:
                    fasta_lengths[current_family_id] = len("".join(current_seq_accumulator).replace("\r", "").replace("\n", "").strip())
                current_family_id = line[1:].strip().split('#')[0].strip()
                current_seq_accumulator = []
            else:
                current_seq_accumulator.append(line)
        if current_family_id and current_seq_accumulator:
            fasta_lengths[current_family_id] = len("".join(current_seq_accumulator).replace("\r", "").replace("\n", "").strip())
            
    merged['Sequence_Length'] = merged['Family_ID'].map(fasta_lengths).fillna(-1).astype(int)
    final_snapshot = merged.fillna('NA')
    
    # Structural Constraint Re-ordering Loop (Locks Column Sledding Grid Anomaly)
    expected_columns = [
        'Family_ID', 'Sequence_Length', 'RepModClass', 'Pfam_Domain', 'Pfam_Evalue', 'DeepTE_Class', 'tRNAscan', 
        'Dfam_Model', 'Dfam_Evalue', 'Dfam_BitScore', 'Heuristic_Target_Protein', 'Fwd_BitScore_Sum', 
        'Fwd_QCovHSP_Max', 'Heuristic_Score', 'RBH_Verdict', 'CD-HIT'
    ]
    df_final = final_snapshot[expected_columns]
    df_final.to_csv(final_curation_tsv, sep='\t', index=False)
    
    print(f"    [✔] Pipeline Matrix Complete. Reference framework saved to '{final_curation_tsv}'.")
    print(df_final['RBH_Verdict'].value_counts().to_string())

# =============================================================================
#                           MAIN INTERFACE CONTROL LOOP
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Unified Upstream Master Pre-Curation Automation Framework (Steps 1-9)")
    parser.add_argument("--query_fna", required=True, help="Path to raw de novo repeat query library FASTA (e.g., mass_scale_2930_te_queries.fna)")
    parser.add_argument("--fwd_hits", required=True, help="Path to raw Diamond forward blastp matches TSV")
    parser.add_argument("--ref_db_fasta", required=True, help="Path to global protein reference database FASTA (e.g., uniref90.fasta)")
    parser.add_argument("--rev_hits", required=True, help="Path to reciprocal tBLASTn reverse mapping hits TSV")
    parser.add_argument("--structural_tsv", required=True, help="Path to structural template baseline matrix TSV")
    parser.add_argument("--cdhit_clstr", required=True, help="Path to CD-HIT sequence clustering log map (.clstr)")
    parser.add_argument("--output_master", default="Hypania_invalida_master_curation.tsv", help="Target filename for compiled final 16-column matrix")
    
    # Internal Pipeline Intermediate Allocations
    parser.add_argument("--tmp_peptides_faa", default="te_viable_peptides.faa", help=argparse.SUPPRESS)
    parser.add_argument("--tmp_fwd_matrix", default="fwd_top10_heuristic_matrix.tsv", help=argparse.SUPPRESS)
    parser.add_argument("--tmp_targets_txt", default="target_uniref_ids.txt", help=argparse.SUPPRESS)
    parser.add_argument("--tmp_reciprocal_faa", default="reciprocal_queries.faa", help=argparse.SUPPRESS)
    parser.add_argument("--tmp_clean_report", default="hypania_cleaned_heuristic_report.tsv", help=argparse.SUPPRESS)
    
    args = parser.parse_args()
    
    print("===============================================================================")
    print("                      INITIALIZING PRE-CURATION PIPELINE                       ")
    print("===============================================================================")
    
    execute_stage1_translation(args.query_fna, args.tmp_peptides_faa, min_orf_len=30)
    
    # NOTE: Downstream alignment calls (DIAMOND/BLAST+) are processed outside this 
    # execution thread via computing clusters to pass output to the following steps.
    
    execute_stage3_heuristic_scoring(args.fwd_hits, args.tmp_fwd_matrix, args.tmp_targets_txt)
    execute_stage4_stream_database(args.ref_db_fasta, args.tmp_targets_txt, args.tmp_reciprocal_faa)
    execute_stage5_6_rbh_and_clean(args.tmp_fwd_matrix, args.rev_hits, args.tmp_clean_report)
    execute_stage7_8_9_consolidation(args.structural_tsv, args.tmp_clean_report, args.cdhit_clstr, args.query_fna, args.output_master)
    
    print("\n[✔] Master Pre-Curation Pipeline Successfully Terminated without Errors.\n")

if __name__ == "__main__":
    main()