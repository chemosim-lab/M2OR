# Utility functions for uniprot database https://www.uniprot.org

from typing import List
from io import StringIO
from urllib import parse
from urllib.request import Request, urlopen
import requests
import pandas

from errors import UniprotNotFoundError, UniprotMultipleOutputError

# def get_uniprot_sequences(uniprot_ids: List, check_consistency: bool = True) -> pandas.DataFrame:
#         """
#         Retrieve uniprot sequences based on a list of uniprot sequence identifier.
# 
#         For large lists it is recommended to perform batch retrieval.
# 
#         Parameters:
#         ----------
#         uniprot_ids: list 
#             list of uniprot identifier
# 
#         Returns:
#         --------
#         df : pandas.DataFrame
#             pandas dataframe with uniprot id column, sequence column and query column.
# 
#         References:
#         -----------
#         Original script:
#         https://www.biostars.org/p/94422/
# 
#         Documentation which columns are available:
#         https://www.uniprot.org/help/uniprotkb%5Fcolumn%5Fnames
#         """
#         url =  'https://www.uniprot.org/uploadlists/'  # This is the webserver to retrieve the Uniprot data
#         params = {
#             'from': "ACC+ID",
#             'to': 'ACC',
#             'format': 'tab',
#             'query': " ".join(uniprot_ids),
#             'columns': 'id,sequence'}
# 
#         data = parse.urlencode(params)
#         data = data.encode('ascii')
#         req = Request(url, data)
#         with urlopen(req) as response:
#             res = response.read()
#         res = res.decode("utf-8")
#         if res == '':
#             raise UniprotNotFoundError('Result from uniprot is empty. Input: {}'.format(uniprot_ids))
#         df_fasta = pandas.read_csv(StringIO(res), sep="\t")
#         df_fasta.columns = ["Entry", "Uniprot_Sequence", "Query"]
#         # it might happen that 2 different ids for a single query id are returned, split these rows
#         df_fasta = df_fasta.assign(Query=df_fasta['Query'].str.split(',')).explode('Query')
#         if check_consistency:
#             set_uniprot_ids = set(uniprot_ids)
#             set_output_ids = set(df_fasta['Entry'])
#             intersect = set_output_ids.intersection(set_uniprot_ids)
#             if len(intersect) < len(set_uniprot_ids):
#                 raise UniprotNotFoundError('Some uniprot IDs were not found: {}'.format(set_uniprot_ids.difference(set_output_ids)))
#             elif len(intersect) > len(set_uniprot_ids):
#                 raise UniprotMultipleOutputError('More uniprot IDs found than inputs. Difference: {}'.format(set_output_ids.difference(set_uniprot_ids)))
#         return df_fasta


def get_uniprot_sequences(uniprot_ids: List, check_consistency: bool = True) -> pandas.DataFrame:
        """
        Retrieve uniprot sequences based on a list of uniprot sequence identifier.

        For large lists it is recommended to perform batch retrieval.

        Parameters:
        ----------
        uniprot_ids: list 
            list of uniprot identifier

        Returns:
        --------
        df : pandas.DataFrame
            pandas dataframe with uniprot id column, sequence column and query column.
        """
        base_url = "https://www.uniprot.org/uniprot/"
        data = []

        for uniprot_id in uniprot_ids:
            response = requests.get(f"{base_url}{uniprot_id}.fasta")
            if response.status_code == 200:
                fasta_data = response.text.split("\n>")
                
                for fasta_entry in fasta_data:
                    if len(fasta_entry) > 0 : 
                        header, sequence = fasta_entry.split("\n", 1)
                    else:
                        print('nothing here')
                    query = header.split("|")[1]
                    data.append({"Entry": uniprot_id, "Uniprot_Sequence": sequence.replace("\n", ""), "Query": query})
            else:
                print(f"Error retrieving sequences, status code: {response.status_code}")
        
        df_fasta = pandas.DataFrame(data)
        # it might happen that 2 different ids for a single query id are returned, split these rows
        df_fasta = df_fasta.assign(Query=df_fasta['Query'].str.split(',')).explode('Query')
        if check_consistency:
            set_uniprot_ids = set(uniprot_ids)
            set_output_ids = set(df_fasta['Entry'])
            intersect = set_output_ids.intersection(set_uniprot_ids)
            if len(intersect) < len(set_uniprot_ids):
                raise UniprotNotFoundError('Some uniprot IDs were not found: {}'.format(set_uniprot_ids.difference(set_output_ids)))
            elif len(intersect) > len(set_uniprot_ids):
                raise UniprotMultipleOutputError('More uniprot IDs found than inputs. Difference: {}'.format(set_output_ids.difference(set_uniprot_ids)))
        return df_fasta


if __name__ == '__main__':
    import os
    df = pandas.read_csv('./RawData/Chemodb_20210223.csv', sep=';', index_col = 0)
    print(df['Uniprot ID'].dropna().unique()[:5])
    df_uniprot = get_uniprot_sequences(df['Uniprot ID'].dropna().unique()[:5])
    df_uniprot = df_uniprot.set_index('Entry', drop = True)
    print(df_uniprot)
    