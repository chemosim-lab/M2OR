import subprocess
import pandas as pd
import numpy as np
import concurrent.futures
from Bio import SeqIO
from typing import List


def blast_search(seq, database_path, blast_executable, species):
    blast_cmd = [blast_executable, "-db", database_path, "-query", "-", "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"]
    blast_output = subprocess.check_output(blast_cmd, input=seq.encode())
    hits = blast_output.decode().strip().split('\n')

    if not hits:
        return float('nan'), float('nan')
    
    hit_data = pd.DataFrame([hit.split('\t') for hit in hits])
    max_evalue = hit_data[11].astype(float).nlargest(2)
    top_hits = hit_data[hit_data[11].astype(float).isin(max_evalue)].copy()


    top_hits['sp_priority'] = top_hits[1].apply(lambda x: x.startswith('sp'))
    
    if species == "homo sapiens":
        top_hits = top_hits.sort_values(by=[5,'sp_priority',8], ascending=[True,False,True])
    else:
        top_hits = top_hits.sort_values(by=[11,5,'sp_priority'], ascending=[False,False,False])
        
    selected_hit = top_hits.iloc[0].copy()
    fasta_id = selected_hit[1]
    uniprotid = selected_hit[1].split('|')[1] if '|' in selected_hit[1] else (selected_hit[1], np.float('nan'))
    identity = float(selected_hit[2])
    return  uniprotid, identity, fasta_id

def search_id_in_fasta(fasta_file, target_id):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return str(fasta_dict[target_id].seq) if target_id in fasta_dict else None


def get_blast_data(mutated_Sequence: List, database_path: List, blast_executable: str, species: List) -> pd.DataFrame:
        """
        Retrieve blast information based on a list of mutated_Sequence.

        
        Parameters:
        ----------
        mutated_Sequence: list 
            list of sequence mutated

        Returns:
        --------
        df : pandas.DataFrame
            pandas dataframe with mutated_Sequence column, blast uniprot id column, blast identity to query sequence column and blast sequence.
        """
        #Subset sequences
        unique_mutated_sequences = mutated_Sequence

        #Parallel blast_search
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(blast_search, unique_mutated_sequences, database_path,
                                        [blast_executable] * len(unique_mutated_sequences), species))

        #Get Results
        df_results = pd.DataFrame(results, columns=["blast_uniprot_id", "blast_identity","blast_fasta_id"])
        df_results["mutated_Sequence"] = unique_mutated_sequences
        df_results["species"] = species


        #Get Sequences
        fasta_file = [f"{x}{'.fasta'}" for x in database_path]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            df_results['blast_seq'] = list(executor.map(search_id_in_fasta, fasta_file, df_results.blast_fasta_id))

        return df_results