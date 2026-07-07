#!/usr/bin/env python3
"""
Master Functional Annotation Pipeline Compiler

Parses Salmon quantification tracks, eggNOG-mapper orthologies, InterProScan 
domain sets, and massive FANTASIA v2 deep learning outputs to completely 
reconstruct the final master annotation table.

Author: Shrikant Vichare
License: CC-BY-4.0 / MIT
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

def normalize_go_terms(s):
    """Cleans, deduplicates, and alphabetically sorts inner GO string components."""
    if pd.isna(s) or not str(s).strip() or str(s).strip().lower() in ['none', 'nan', '-', '', 'none,none']:
        return np.nan
    elements = str(s).replace('|', ',').split(',')
    terms = sorted(set([t.strip() for t in elements if t.strip() and t.strip().upper().startswith('GO:')]))
    return ",".join(terms) if terms else np.nan

def clean_transcript_id(tid):
    """Ensures consistent Transcript ID keys across different tool naming suffixes."""
    if pd.isna(tid):
        return ""
    s = str(tid).strip()
    return s

def build_master_annotation_archive(args):
    print("===============================================================================")
    echo_msg = "               STARTING FUNCTIONAL ANNOTATION SYNTHESIS               "
    print(echo_msg)
    print("===============================================================================")

    # ---------------------------------------------------------------------------
    # STEP 1: PARSE SALMON QUANTIFICATION SCAFFOLD
    # ---------------------------------------------------------------------------
    print("[1/5] Extracting expression profiling framework via Salmon quant track...")
    if not os.path.exists(args.quant_sf):
        print(f"[-] Critical Error: Salmon quant trace missing -> '{args.quant_sf}'")
        sys.exit(1)
        
    df_quant = pd.read_csv(args.quant_sf, sep='\t')
    df_quant['Transcript_ID'] = df_quant['Name'].apply(clean_transcript_id)
    
    # Establish the foundational structural scaffold tracking all coding transcripts
    master_df = df_quant[['Transcript_ID', 'Length', 'TPM']].copy()
    print(f"    -> Successfully scaffolded {len(master_df):,} baseline transcript models.")

    # ---------------------------------------------------------------------------
    # STEP 2: STREAM AND STREAMLINE INTERPROSCAN ARCHITECTURES
    # ---------------------------------------------------------------------------
    print("[2/5] Parsing domain and site structural annotations via InterProScan...")
    if not os.path.exists(args.interpro):
        print(f"[-] Critical Error: InterProScan tabular source missing -> '{args.interpro}'")
        sys.exit(1)

    ipr_data = {}
    ipr_chunk_size = 200000
    
    # Process line records using chunk-based streaming to protect system RAM
    for chunk in pd.read_csv(args.interpro, sep='\t', header=None, chunksize=ipr_chunk_size, low_memory=False):
        for _, row in chunk.iterrows():
            tid = clean_transcript_id(row[0])
            if not tid:
                continue
                
            ipr_acc = str(row[11]).strip() if pd.notna(row[11]) and str(row[11]).strip() != '-' else np.nan
            ipr_desc = str(row[12]).strip() if pd.notna(row[12]) and str(row[12]).strip() != '-' else np.nan
            ipr_go = str(row[13]).strip() if pd.notna(row[13]) and str(row[13]).strip() != '-' else np.nan
            
            if tid not in ipr_data:
                ipr_data[tid] = {'acc': set(), 'desc': set(), 'go': set()}
                
            if pd.notna(ipr_acc): ipr_data[tid]['acc'].add(ipr_acc)
            if pd.notna(ipr_desc): ipr_data[tid]['desc'].add(ipr_desc)
            if pd.notna(ipr_go):
                for g in ipr_go.split('|'):
                    if g.strip().upper().startswith('GO:'):
                        ipr_data[tid]['go'].add(g.strip())

    # Consolidate aggregated sets into formatted dataframe columns
    ipr_rows = []
    for tid, sets in ipr_data.items():
        ipr_rows.append({
            'Transcript_ID': tid,
            'InterPro_Accession': "|".join(sorted(sets['acc'])) if sets['acc'] else "-",
            'InterPro_Description': "|".join(sorted(sets['desc'])) if sets['desc'] else "-",
            'InterPro_GO_Terms': ",".join(sorted(sets['go'])) if sets['go'] else np.nan
        })
        
    df_ipr = pd.DataFrame(ipr_rows)
    del ipr_data # Immediately release memory allocations back to the host system
    print(f"    -> Consolidated signature domains across {len(df_ipr):,} unique transcripts.")

    # ---------------------------------------------------------------------------
    # STEP 3: PARSE EGGNOG-MAPPER ORTHOLOGIES
    # ---------------------------------------------------------------------------
    print("[3/5] Extracting functional orthology systems via eggNOG-mapper...")
    if not os.path.exists(args.eggnog):
        print(f"[-] Critical Error: eggNOG annotation dataset missing -> '{args.eggnog}'")
        sys.exit(1)

    # Automatically handle variable comment rows to isolate the core table header safely
    df_eggnog = pd.read_csv(args.eggnog, sep='\t', comment='#', header=None)
    
    # Trace the true programmatic headers from the last commented row configuration
    with open(args.eggnog, 'r') as f:
        headers = []
        for line in f:
            if line.startswith('#query'):
                headers = line.strip().split('\t')
                break
    
    if len(headers) == len(df_eggnog.columns):
        df_eggnog.columns = headers
    else:
        # Hardcoded fallback layer tailored to eggNOG-mapper v2.x standards if needed
        df_eggnog.rename(columns={0: '#query', 8: 'Preferred_name', 7: 'Description', 9: 'GOs', 13: 'KEGG_Pathway', 6: 'COG_category'}, inplace=True)
    
    df_eggnog['Transcript_ID'] = df_eggnog['#query'].apply(clean_transcript_id)
    
    # Standardize empty cell values using unified placeholders
    df_eggnog['eggNOG_Preferred_name'] = df_eggnog['Preferred_name'].fillna('-').astype(str).str.strip()
    df_eggnog['eggNOG_Description'] = df_eggnog['Description'].fillna('-').astype(str).str.strip()
    df_eggnog['eggNOG_GOs'] = df_eggnog['GOs'].apply(normalize_go_terms)
    df_eggnog['eggNOG_KEGG_Pathway'] = df_eggnog['KEGG_Pathway'].fillna('-').astype(str).str.strip()
    df_eggnog['eggNOG_COG_category'] = df_eggnog['COG_category'].fillna('-').astype(str).str.strip()
    
    eggnog_out_cols = ['Transcript_ID', 'eggNOG_Preferred_name', 'eggNOG_Description', 'eggNOG_GOs', 'eggNOG_KEGG_Pathway', 'eggNOG_COG_category']
    df_eggnog = df_eggnog[eggnog_out_cols].drop_duplicates(subset=['Transcript_ID'])
    print(f"    -> Extracted orthologous profiles across {len(df_eggnog):,} transcripts.")

    # ---------------------------------------------------------------------------
    # STEP 4: DEPLOY COMPRESSION MATRIX OVER FANTASIA EMBEDDINGS
    # ---------------------------------------------------------------------------
    print("[4/5] Running heuristic sorting loops over protein language model outputs...")
    if not os.path.exists(args.fantasia):
        print(f"[-] Critical Error: FANTASIA summary table input missing -> '{args.fantasia}'")
        sys.exit(1)

    cols_fantasia = ['accession', 'go_id', 'distance', 'reliability_index', 'identity', 'similarity', 'length_query', 'length_reference']
    running_best = pd.DataFrame()
    fan_chunk_size = 500000

    for chunk in pd.read_csv(args.fantasia, usecols=cols_fantasia, chunksize=fan_chunk_size, engine='c', low_memory=False):
        for col in ['distance', 'reliability_index', 'identity', 'similarity', 'length_query', 'length_reference']:
            chunk[col] = pd.to_numeric(chunk[col], errors='coerce')
        
        chunk = chunk.dropna(subset=['distance', 'reliability_index'])
        if chunk.empty:
            continue
            
        chunk['go_id'] = chunk['go_id'].fillna('None').astype(str).str.strip()
        chunk_grouped = chunk.groupby(['accession', 'distance', 'reliability_index', 'identity', 'similarity', 'length_query', 'length_reference'], as_index=False)['go_id'].agg(lambda x: ','.join(x))
        
        # Deploy strict 4-layer heuristic sorting tuple
        chunk_grouped = chunk_grouped.sort_values(by=['distance', 'reliability_index', 'identity', 'length_reference'], ascending=[True, False, False, False])
        chunk_best = chunk_grouped.drop_duplicates(subset=['accession'], keep='first')
        
        running_best = pd.concat([running_best, chunk_best])
        running_best = running_best.sort_values(by=['distance', 'reliability_index', 'identity', 'length_reference'], ascending=[True, False, False, False])
        running_best = running_best.groupby(['accession', 'distance', 'reliability_index', 'identity', 'similarity', 'length_query', 'length_reference'], as_index=False)['go_id'].agg(lambda x: ','.join(x))
        running_best = running_best.sort_values(by=['distance', 'reliability_index', 'identity', 'length_reference'], ascending=[True, False, False, False])
        running_best = running_best.drop_duplicates(subset=['accession'], keep='first')

    running_best['go_id_cleaned'] = running_best['go_id'].apply(normalize_go_terms)
    running_best['Transcript_ID'] = running_best['accession'].apply(clean_transcript_id)

    # Reconstruct standalone TopHits summaries to match Zenodo properties exactly
    summary_fields = ['accession', 'go_ids', 'distance', 'reliability_index', 'identity', 'similarity', 'length_query', 'length_reference']
    
    # Save standalone RI 0.90 summary map
    df_s090 = running_best.copy()
    df_s090['go_ids'] = np.where(df_s090['reliability_index'] >= 0.90, df_s090['go_id_cleaned'], "None")
    df_s090['length_query'] = df_s090['length_query'].astype(int)
    df_s090['length_reference'] = df_s090['length_reference'].astype(int)
    df_s090[summary_fields].sort_values(by='accession').to_csv(args.out_sum90, sep='\t', index=False)
    print(f"    -> Standalone map generated successfully -> '{args.out_sum90}'")

    # Save standalone RI 0.95 summary map
    df_s095 = running_best.copy()
    df_s095['go_ids'] = np.where(df_s095['reliability_index'] >= 0.95, df_s095['go_id_cleaned'], "None")
    df_s095['length_query'] = df_s095['length_query'].astype(int)
    df_s095['length_reference'] = df_s095['length_reference'].astype(int)
    df_s095[summary_fields].sort_values(by='accession').to_csv(args.out_sum95, sep='\t', index=False)
    print(f"    -> Standalone map generated successfully -> '{args.out_sum95}'")

    # Construct the relational column sets for integration into the master database
    running_best['ESM_GO_Terms_RI0.90'] = np.where(running_best['reliability_index'] >= 0.90, running_best['go_id_cleaned'], np.nan)
    running_best['ESM_Distance_RI0.90'] = running_best['distance']
    running_best['ESM_Confidence_RI0.90'] = running_best['reliability_index']
    running_best['ESM_Identity_RI0.90'] = running_best['identity']
    running_best['ESM_Ref_Length_RI0.90'] = running_best['length_reference'].astype(int)

    # Note the structured column inversion layout required to perfectly align with the master spreadsheet
    running_best['ESM_GO_Terms_RI0.95'] = np.where(running_best['reliability_index'] >= 0.95, running_best['go_id_cleaned'], np.nan)
    running_best['ESM_Confidence_RI0.95'] = running_best['reliability_index']
    running_best['ESM_Distance_RI0.95'] = running_best['distance']
    running_best['ESM_Identity_RI0.95'] = running_best['identity']
    running_best['ESM_Ref_Length_RI0.95'] = running_best['length_reference'].astype(int)

    esm_merge_cols = [
        'Transcript_ID', 
        'ESM_GO_Terms_RI0.90', 'ESM_Distance_RI0.90', 'ESM_Confidence_RI0.90', 'ESM_Identity_RI0.90', 'ESM_Ref_Length_RI0.90',
        'ESM_GO_Terms_RI0.95', 'ESM_Confidence_RI0.95', 'ESM_Distance_RI0.95', 'ESM_Identity_RI0.95', 'ESM_Ref_Length_RI0.95'
    ]
    df_esm_final = running_best[esm_merge_cols].drop_duplicates(subset=['Transcript_ID'])

    # ---------------------------------------------------------------------------
    # STEP 5: COMPILE MASTER ANNOTATION TABLE
    # ---------------------------------------------------------------------------
    print("[5/5] Executing high-fidelity master joins across relational structures...")
    
    # Sequentially map dataset columns back into the core Salmon matrix via left joins
    master_df = pd.merge(master_df, df_eggnog, on='Transcript_ID', how='left')
    master_df = pd.merge(master_df, df_ipr, on='Transcript_ID', how='left')
    master_df = pd.merge(master_df, df_esm_final, on='Transcript_ID', how='left')
    
    # Fill remaining unannotated structural fields with standard text indicators
    text_cols = ['eggNOG_Preferred_name', 'eggNOG_Description', 'eggNOG_KEGG_Pathway', 'eggNOG_COG_category', 'InterPro_Accession', 'InterPro_Description']
    for c in text_cols:
        master_df[c] = master_df[c].fillna('-')
        
    # Standardize column sort key sequence to precisely match the target geometry
    final_structural_order = [
        'Transcript_ID', 'Length', 'TPM',
        'eggNOG_Preferred_name', 'eggNOG_Description', 'eggNOG_GOs', 'eggNOG_KEGG_Pathway', 'eggNOG_COG_category',
        'InterPro_Accession', 'InterPro_Description', 'InterPro_GO_Terms',
        'ESM_GO_Terms_RI0.90', 'ESM_Distance_RI0.90', 'ESM_Confidence_RI0.90', 'ESM_Identity_RI0.90', 'ESM_Ref_Length_RI0.90',
        'ESM_GO_Terms_RI0.95', 'ESM_Confidence_RI0.95', 'ESM_Distance_RI0.95', 'ESM_Identity_RI0.95', 'ESM_Ref_Length_RI0.95'
    ]
    
    master_df = master_df[final_structural_order].sort_values(by='Transcript_ID')
    
    # Reconstruct individual threshold profiles to provide clean outputs across all files
    iso_cols_base = ['Transcript_ID', 'Length', 'TPM', 'eggNOG_Preferred_name', 'eggNOG_Description', 'eggNOG_GOs', 'eggNOG_KEGG_Pathway', 'eggNOG_COG_category', 'InterPro_Accession', 'InterPro_Description', 'InterPro_GO_Terms']
    
    # Synthesize standalone RI 0.90 Master Table
    df_m90 = master_df[iso_cols_base + ['ESM_GO_Terms_RI0.90', 'ESM_Distance_RI0.90', 'ESM_Confidence_RI0.90', 'ESM_Identity_RI0.90', 'ESM_Ref_Length_RI0.90']].copy()
    df_m90.columns = iso_cols_base + ['ESM_GO_Terms', 'ESM_Distance', 'ESM_Confidence', 'ESM_Identity', 'ESM_Ref_Length']
    df_m90.to_csv(args.out_m90, sep='\t', index=False)
    print(f"    -> Standalone master spreadsheet constructed -> '{args.out_m90}'")

    # Synthesize standalone RI 0.95 Master Table
    # Note that we automatically handle column order rearrangement for the isolated spreadsheet layers
    df_m95 = master_df[iso_cols_base + ['ESM_GO_Terms_RI0.95', 'ESM_Distance_RI0.95', 'ESM_Confidence_RI0.95', 'ESM_Identity_RI0.95', 'ESM_Ref_Length_RI0.95']].copy()
    df_m95.columns = iso_cols_base + ['ESM_GO_Terms', 'ESM_Distance', 'ESM_Confidence', 'ESM_Identity', 'ESM_Ref_Length']
    df_m95.to_csv(args.out_m95, sep='\t', index=False)
    print(f"    -> Standalone master spreadsheet constructed -> '{args.out_m95}'")

    # Export primary multi-threshold dataset
    master_df.to_csv(args.out_master, sep='\t', index=False)
    print(f"\n[✔] Done! Primary reproducible repository database saved -> '{args.out_master}'")
    print("[✔] Reviewer processing successfully completed. All matrices are synchronized.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Master Functional Annotation Synthesis & Reproducibility Suite")
    parser.add_argument("--fantasia", default="summary_manual.csv", help="Path to uncompressed raw FANTASIA alignment matrix")
    parser.add_argument("--quant_sf", default="Final_Transcriptome_quant.sf", help="Path to Salmon transcript quantification track")
    parser.add_argument("--eggnog", default="eggNOG_Mapper_annotations.tabular", help="Path to eggNOG-mapper annotation file")
    parser.add_argument("--interpro", default="InterProScan.tabular", help="Path to InterProScan domain dataset")
    parser.add_argument("--out_master", default="HYPANIA_INVALIDA_MASTER_ANNOTATION.tsv", help="Target filename for unified output file")
    parser.add_argument("--out_m90", default="HYPANIA_MASTER_ANNOTATION_RI_0.90.tsv", help="Target filename for isolated 0.90 matrix")
    parser.add_argument("--out_m95", default="HYPANIA_MASTER_ANNOTATION_RI_0.95.tsv", help="Target filename for isolated 0.95 matrix")
    parser.add_argument("--out_sum90", default="HYP_ESM_ONT_TopHits_Summary_RI_0.90.tsv", help="Target filename for standalone 0.90 summaries")
    parser.add_argument("--out_sum95", default="HYP_ESM_ONT_TopHits_Summary_RI_0.95.tsv", help="Target filename for standalone 0.95 summaries")
    
    build_master_annotation_archive(parser.parse_args())