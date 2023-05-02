from errors import MutationError
import pandas

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def _perform_mutation(mutations, seq):
    """
    Main logic for performing mutations.
    """    
    for mutation in mutations:
        # print('->' + mutation)
        _from = mutation[0]
        _to = mutation[-1]
        _position = int(mutation[1:-1]) - 1
        if seq[_position] == _from:
            seq = seq[:_position] + _to + seq[_position+1:]
        else:
            left_pos = _position - 5
            right_pos = _position + 4
            if _position < 5: left_pos = 0
            if _position > len(seq) - 4: right_pos = len(seq)
            raise MutationError('Expected letter {} on position {} in sequence arround: {}. Found {}. Mutation: {}'.format(_from, _position, seq[left_pos:right_pos], seq[_position], mutation))
    return seq

def perform_mutation(row, mutation_col = 'Mutation', seq_col = 'Sequence'):
    """
    Perform mutation on a row in pandas.DataFrame. This function is designed to be used in pandas apply.

    Paramters:
    ----------
    row : pandas.Series
        row of pandas dataframe

    mutation_col : str
        name of the column with mutation information.

    seq_col : str
        name of the column with sequences.

    Returns:
    --------
    mutated_seq : str
        mutated sequence
    """
    try:
        if row[mutation_col] == row[mutation_col] and row[seq_col] == row[seq_col]:
            seq = row[seq_col]
            mutations = row[mutation_col].strip().split('_')
            return _perform_mutation(mutations, seq)
        else:
            return row[seq_col]
    except MutationError as e:
        # print(e)
        raise e
        # return float('nan')
    except Exception as e:
        print(row)
        print(row[seq_col])
        raise  e


def merge_cols_with_priority(row, primary_col = 'Uniprot_Sequence', secondary_col = 'Sequence'): # merge_seqs
    """
    merge two columns. primary column is kept and only if entry is missing use secondary_col.
    This should be used in pandas apply.

    Paramters:
    ----------
    row : pandas.Series
        row of pandas dataframe.
    
    primary_col : str
        name of the main column.

    secondary_col : str
        name of the secondary column used only when info in the first is missing.

    Returns:
    --------
    entry : Any
        value to put to the merged column.
    """
    if row[primary_col] == row[primary_col]:
        return row[primary_col]
    else:
        return row[secondary_col]


def enumerate_isomers(canonicalSMILES):
    """
    Get ismoers for a given canonical SMILES.

    Parameters:
    -----------
    canonicalSMILES : str
        canonical SMILES.

    Returns:
    --------
    isomers : list
        list of isomers found by rdkit EnumerateStereoisomers.

    References:
    -----------
    https://www.rdkit.org/docs/source/rdkit.Chem.EnumerateStereoisomers.html
    """
    mol = Chem.MolFromSmiles(canonicalSMILES)
    isomers = tuple(EnumerateStereoisomers(mol))
    isomers = [Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers]
    return isomers


if __name__ == '__main__':
    print(enumerate_isomers('CCC(C)CO.CC(C)CCO'))
