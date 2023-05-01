from unittest import result
import pandas
import pubchempy
import json
from urllib.error import HTTPError
import time

def get_map_inchikey_to_CID(inchikey):
    mols = pubchempy.get_compounds(inchikey, 'inchikey')
    return {mol.inchikey: mol.cid for mol in mols if mol is not None}


def get_map_inchikey_to_synonyms(inchikey):
    mols = pubchempy.get_compounds(inchikey, 'inchikey')
    result = {}
    for mol in mols:
        if mol is not None:
            if mol.iupac_name is not None:
                result[mol.inchikey] = mol.synonyms + [mol.iupac_name]
            else:
                result[mol.inchikey] = mol.synonyms
    return result


def get_map_inchikey_to_canonicalSMILES(inchikey):
    mols = pubchempy.get_compounds(inchikey, namespace = u'inchikey')
    return {mol.inchikey: mol.canonical_smiles for mol in mols if mol is not None}


def get_map_isomericSMILES_to_inchikey(smiles):
    if isinstance(smiles, str):
        smiles = [smiles]
    res = {}
    i = 0
    while i < len(smiles):
        try:
            response = pubchempy.get_compounds(smiles[i], 'smiles')
            if len(response) == 1:
                mol = response[0]
                res[smiles[i]] = mol.inchikey
            elif len(response) > 1:
                raise ValueError('WARNING: More than one InChI Key for smiles: {}'.format(smiles[i]))
            else:
                raise ValueError('WARNING: No InChIKey found')
        except HTTPError:
            time.sleep(1)
        i += 1
    return res


def get_map_name_to_inchikeys(name):
    if isinstance(name, str):
        name = [name]
    res = {}
    i = 0
    while i < len(name):
        try:
            response = pubchempy.get_compounds(name[i], 'name')
            if len(response) >= 1:
                res[name[i]] = [mol.inchikey for mol in response]
        except HTTPError:
            time.sleep(1)
        i += 1
    return res




if __name__ == '__main__':
    import os
    df = pandas.read_csv('/home/matej/Desktop/tmp/Chemodb_20210110.csv', sep=';', index_col = 0)
    inchikey = df['InchI Key'].str.strip().str.replace('\s\s+', ' ').str.split(' ').explode() # TODO: pandas FutureWarning for this row.
    map_inchikey_to_CID = get_map_inchikey_to_CID(list(inchikey.unique()))

    df_map_inchikey_to_CID = pandas.Series(map_inchikey_to_CID)
    df_map_inchikey_to_CID.index.name = 'InChI Key'
    df_map_inchikey_to_CID.name = 'CID'
    df_map_inchikey_to_CID.to_csv(os.path.join('Data', 'map_inchikey_to_CID.csv'), sep = ';')
    
