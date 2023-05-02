# Formating before merging and checking.
import os
import sys
import pandas
import re
import json
import itertools
import logging


from utils import perform_mutation, merge_cols_with_priority, enumerate_isomers
from uniprot_utils import get_uniprot_sequences
from blast_utils import get_blast_data
from pubchem_utils import get_map_inchikey_to_canonicalSMILES, get_map_isomericSMILES_to_inchikey

# (OK) TODO: Order mutations (for mutated_Uniprot_ID)
# (OK) TODO: Stip spaces
# (OK) TODO: Lowercase Mixture name, type, ...
# (OK) TODO: luciferase assay vs luciferase
# (OK) TODO: normalize separator
# (OK) TODO: exclude examples without sequence to separate csv.


class PreFormatter:
    """
    Perform basic formatting like white space normalizations and deletions, lowercasing and replacing obvious mistakes.

    This class also excludes rows without sequence or with exotic classes like antagonists to a separate dataframe.

    Attributes:
    -----------
    auxillary_dir : str
        directory with auxillary data like \'uniprot_sequences.csv\'.

    log_dir : str
        loggig dir name.

    logger : logging.Logger
        logger
    
    df_uniprot_cols : list
        columns expected to be found in \'uniprot_sequences.csv\'. This serves as a precaution.
    
    strip_whitespace_cols : list
        list of columns to which we apply strip_whitespace method.

    lowercase_string_cols : list
        list of columns to which we apply lowercase_string method.
    
    replace_cols : list
        list of dictionaries in the form {'col' : colname, 'from' : list_of_patterns, 'to' : pattern}. This is 
        used in replace_ method.
    
    normalize_entries_col : list
        list of dictionaries in the form {'col' : colname, 'from' : words_to_replace, 'to' : replace_with_this}. This is 
        used in normalize_entries method to fix obvious problems like writing \'luciferase assay\' instead of \'luciferase\'.

    conditioned_set_value_cols : list
        list of dictionaries in the form {'cond_col' : colname, 'cond_val' : value, 'target_col' : colname, 'target_val' : value} which is
        used to set values for column \'target_col\' in entries where \'cond_col\' are \'cond_val\' to \'target_val\'. For instance if there is 
        \'n.d\' in Value column we set Unit to \'n.d\'

    df_uniprot : pandas.DataFrame
        auxillary dataframe with mapping from uniprot ID to sequence.
    """
    def __init__(self, auxillary_dir = 'Data', log_dir = 'logs'):
        self.auxillary_dir = auxillary_dir
        self.log_dir = log_dir

        if not os.path.exists(self.auxillary_dir):
            os.makedirs(self.auxillary_dir)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        self.logger = logging.getLogger(__class__.__name__)
        self.logger.setLevel(logging.DEBUG)
        
        logger_formatter = logging.Formatter("%(levelname)-10s:\t%(asctime)-20s:\t%(message)s", "%Y-%m-%d %H:%M:%S")
        logger_file_handler = logging.FileHandler(os.path.join(self.log_dir, __class__.__name__ + '.log'), mode = 'w')
        logger_file_handler.setFormatter(logger_formatter)
        logger_stdout_handler = logging.StreamHandler(sys.stdout)

        logger_file_handler.setLevel(logging.DEBUG)
        logger_stdout_handler.setLevel(logging.INFO)

        self.logger.addHandler(logger_file_handler)
        self.logger.addHandler(logger_stdout_handler)

        self.df_uniprot_cols = ["Entry", "Uniprot_Sequence", "Query"]

        self.strip_whitespace_cols = ['species', 'Mutation', 'Gene ID', 'Uniprot ID', 'Sequence', 'Name',
                                    'CID', 'CAS', 'InChI Key', 'canonicalSMILES', 'Parameter', 'Value', 'Unit', 'Value_Screen',
                                    'Unit_Screen', 'Responsive', 'nbr_measurements', 'Type', 'Cell_line', 'Co_transfection', 'Assay System', 
                                    'Tag', 'Reference', 'DOI', 'Reference Position', 'Mixture']

        self.lowercase_string_cols = ['species', 'Parameter', 'Value', 'Value_Screen', 'Co_transfection', 'Assay System','Gprotein', 'Delivery', 'Assay', 'Tag', 'Mixture']

        self.replace_cols = [{'col' : 'Mutation', 'from' : [' ', '-'], 'to' : '_'},
                            {'col' : 'InChI Key', 'from' : ['_'], 'to' : ' '},
                            {'col' : 'Sequence', 'from' : ['\n', ' '], 'to' : ''},
                            {'col' : 'Value', 'from' : ['> '], 'to' : '>'},
                            {'col' : 'Value', 'from' : ['< '], 'to' : '<'},
                            ]

        self.conditioned_set_value_cols = [{'cond_col' : 'Value', 'cond_val' : 'n.d', 'target_col' : 'Unit', 'target_val' : 'n.d'},
                                        ]

        self.normalize_entries_col = [{'col' : 'Type', 'from' : 'luciferase assay', 'to' : 'luciferase'}, # TODO maybe change this ??
                                    {'col' : 'Value', 'from' : 'n.d.', 'to' : 'n.d'},
        ]


    def _load_auxillary(self):
        """
        load auxillary files: \'uniprot_sequences.csv\' and put the result to attributes.

        Returns:
        --------
        df_uniprot : pandas.DataFrame
            loaded \'uniprot_sequences.csv\'.
        """
        try:
            df_uniprot = pandas.read_csv(os.path.join(self.auxillary_dir, 'uniprot_sequences.csv'), sep = ';', index_col = None)
        except FileNotFoundError:
            df_uniprot = pandas.DataFrame([], columns = self.df_uniprot_cols)
            # raise FileNotFoundError('Support file: \'uniprot_sequences.csv\' not found. Please call self.update_df_uniprot or run uniprot_utils.py to download sequences.')

        assert len(df_uniprot.columns) == len(self.df_uniprot_cols)
        for i in range(len(df_uniprot.columns)):
            if df_uniprot.columns[i] != self.df_uniprot_cols[i]:
                raise ValueError('df_uniprot has different columns or column positions than in self.df_uniprot_cols')
        df_uniprot.set_index(self.df_uniprot_cols[0], drop = True, inplace = True)
        self.df_uniprot = df_uniprot
        return df_uniprot  


    def _update_auxillary_df_uniprot(self, full_df):
        """
        update \'uniprot_sequences.csv\' with new uniprot sequences found in full_df.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'Uniprot ID\' columns which will be used to check which IDs are new.
        
        Returns:
        --------
        df_uniprot : padnas.DataFrame
            updated df_uniprot
        """
        df_uniprot = self.df_uniprot.copy()
        candidate_idx = pandas.Index(full_df['Uniprot ID'].dropna().unique())
        new_idx = candidate_idx.difference(df_uniprot.index)
        if len(new_idx) > 0:
            self.logger.info('Updating df_uniprot...')
            NEW = get_uniprot_sequences(new_idx.tolist())
            NEW.set_index(self.df_uniprot_cols[0], drop = True, inplace = True)
            df_uniprot = df_uniprot.append(NEW, ignore_index = False, verify_integrity = True)
            df_uniprot.to_csv(os.path.join(self.auxillary_dir, 'uniprot_sequences.csv'), sep = ';', index = True)
        return df_uniprot

    def _update_auxillary(self, full_df):
        """
        call to all _update_auxillary_*.
        """
        self.df_uniprot = self._update_auxillary_df_uniprot(full_df)
        return 


    def add_implied_columns(self, df):
        """
        Add \'_Sequence\' column to df with retrieved sequences from UniProt (or keeping sequence if no UniProt ID is provided).

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        df_join : pandas.DataFrame
            dataframe with added column.
        """
        self._load_auxillary()
        self._update_auxillary(df)
        df_join = df.copy()
        if df_join["Uniprot ID"].isna().all():
            df_join['Uniprot_Sequence'] = float('nan')
            #print('empty')
        else:
            print(df_join.columns)
            df_join = df_join.join(self.df_uniprot[['Uniprot_Sequence']], on = 'Uniprot ID', how = 'left')
        
        df_join['_Sequence'] = df_join.apply(lambda x: merge_cols_with_priority(x, primary_col = 'Uniprot_Sequence', secondary_col = 'Sequence'), axis = 1)
        return df_join
        
    @staticmethod
    def _strip_whitespace(x):
        if isinstance(x, str):
            return re.sub(r'\s\s+',' ', x.strip()) # .replace('\s\s+', ' ').replace(' +', ' ')
        else:
            return x

    def strip_whitespace(self, df):
        """
        for each entry in \'strip_whitespace_cols\' normalize whitespaces.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with updated columns.
        """
        transformed_df = df.copy()
        for col in self.strip_whitespace_cols:
            transformed_df[col] = transformed_df[col].apply(self._strip_whitespace)
        return transformed_df


    @staticmethod
    def _lower(x):
        if isinstance(x, str):
            return x.lower()
        else:
            return x

    def lowercase_string(self, df):
        """
        for each entry in \'lowercase_string_cols\' put entries to lowercase.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with updated columns.
        """
        transformed_df = df.copy()
        for col in self.lowercase_string_cols:
            transformed_df[col] = transformed_df[col].apply(self._lower)
        return transformed_df


    @staticmethod
    def _replace_(x, _from, _to):
        if isinstance(x, str):
            for s in _from:
                x = x.replace(s, _to)
            return x
        else:
            return x

    def replace_(self, df):
        """
        for each entry in \'replace_cols\', for a for a column given by \'col\' replace patterns in \'from\' to pattern in \'to\'.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with updated columns.
        """
        transformed_df = df.copy()
        for case in self.replace_cols:
            transformed_df[case['col']] = transformed_df[case['col']].apply(lambda x: self._replace_(x, case['from'], case['to']))
        return transformed_df


    @staticmethod
    def _normalize_entries(x, from_entry, to_entry):
        if isinstance(x, str):
            return x.replace(from_entry, to_entry)
        else:
            return x

    def normalize_entries(self, df):
        """
        for each entry in \'normalize_entries_col\', for a column given by \'col\' replace phrase in \'from\' to phrase in \'to\'.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with updated columns.
        """
        transformed_df = df.copy()
        for case in self.normalize_entries_col:
            transformed_df[case['col']] = transformed_df[case['col']].apply(lambda x: self._normalize_entries(x, case['from'], case['to']))
        return transformed_df

    
    def conditioned_set_value(self, df):
        """
        For each entry in \'conditioned_set_value_cols\', for entries where value in \'cond_col\' is \'cond_val\' set value of another 
        column \'target_col\' to \'target_val\'.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with updated columns.
        """
        transformed_df = df.copy()
        for case in self.conditioned_set_value_cols:
            _idx = transformed_df[transformed_df[case['cond_col']] == case['cond_val']].index
            transformed_df.loc[_idx, case['target_col']] = case['target_val']
        return transformed_df


    def exclude_empty_sequence(self, df, excluded_df = None):
        """
        Exclude records without sequence to a separate dataframe. 
        Number of excluded rows and current size of exclude_df is logged.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        excluded_df : pandas.DataFrame, optional (default = None)
            previously excluded entries. Newly excluded are added to this dataframe.
            If None, an empty dataframe with the same structure as df is created.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with all the entries having sequence available.
        
        excluded_df : pandas.DataFrame
            remaining of df. We are not able to retrieve sequence for these entries.
        """
        transformed_df = df.copy()
        if excluded_df is None:
            excluded_df = df.iloc[:0,:].copy() # Empty dataframe with same structure as df.
        new_excluded_df = excluded_df.copy()
        df_join = self.add_implied_columns(df)
        discard_idx = df_join[df_join['_Sequence'].isna()].index
        
        # excluded_df = excluded_df.loc[discard_idx]
        # transformed_df = transformed_df.loc[transformed_df.index.difference(discard_idx)]

        new_excluded_df = new_excluded_df.append(transformed_df.loc[discard_idx])
        transformed_df = transformed_df.loc[transformed_df.index.difference(discard_idx)] #If excluded df empty create empty dataframe

        self.logger.info('Num excluded: {}'.format(len(discard_idx)))
        self.logger.info('Exclude_df size: {}'.format(len(new_excluded_df)))
        return transformed_df, new_excluded_df

    def exclude_antagonist(self, df, excluded_df = None): 
        '''
        Exclude records with exotic class to excluded_df dataframe.
        Number of excluded rows and current size of exclude_df is logged.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        excluded_df : pandas.DataFrame, optional (default = None)
            previously excluded entries. Newly excluded are added to this dataframe.
            If None, an empty dataframe with the same structure as df is created.

        Returns:
        --------
        transformed_df : pandas.DataFrame
            copy of df with all the entries having sequence available.
        
        new_excluded_df : pandas.DataFrame
            excluded rows appended to excluded_df. These rows have some exotic class.

        Notes:
        ------
        TODO Change the way to exclude antagonist and inverse agonists. Create another column for antagonist ? 
        '''
        transformed_df = df.copy()
        if excluded_df is None:
            excluded_df = df.iloc[:0,:].copy() # Empty dataframe with same structure as df.
        new_excluded_df = excluded_df.copy()
        condition = transformed_df['Responsive'].apply(lambda x: True if re.match("2|-1", str(x)) else False)
        discard_idx = transformed_df[condition].index
        new_excluded_df = new_excluded_df.append(transformed_df.loc[discard_idx])
        transformed_df = transformed_df.loc[transformed_df.index.difference(discard_idx)] #If excluded df empty create empty dataframe
        self.logger.info('Num excluded: {}'.format(len(discard_idx)))
        self.logger.info('Exclude_df size: {}'.format(len(new_excluded_df)))
        return transformed_df, new_excluded_df

    def __call__(self, df):
        """
        Paramters:
        ----------
        df : pandas.DataFrame
            raw dataframe.

        Returns:
        --------
        df : pandas.DataFrame
            formatted dataframe. All records here have sequence available and are either \'agonist\' or \'non-binding\'.

        exclude_df : pandas.DataFrame
        excluded records either without sequence or having some exotic class.
        """
        df = self.strip_whitespace(df)
        df = self.lowercase_string(df)
        df = self.replace_(df)
        df = self.normalize_entries(df)
        df = self.conditioned_set_value(df)
        df, excluded_df = self.exclude_empty_sequence(df)
        df, excluded_df = self.exclude_antagonist(df, excluded_df)
        return df, excluded_df



class PostFormatter:
    """
    Perform formatting after checker. This serves for adding new implied columns such as \'Mixtures\' or \'mutated_Sequence\'.

    Attributes:
    -----------
    auxillary_dir : str
        directory with auxillary data like \'map_inchikey_to_canonicalSMILES.json\' and \'uniprot_sequences.csv\'.

    log_dir : str
        loggig dir name.

    logger : logging.Logger
        logger

    df_uniprot_cols : list
        columns expected to be found in \'uniprot_sequences.csv\'. This serves as a precaution.

    map_inchikey_to_canonicalSMILES : dict
        auxillary dictionary mapping InChI key to canonical SMILES.
    """
    def __init__(self, log_to_file = True, auxillary_dir = None, log_dir = 'logs'):
        if auxillary_dir is None:
            self.auxillary_dir = 'Data'
        else:
            self.auxillary_dir = auxillary_dir
        self.log_dir = log_dir

        self.logger = logging.getLogger(__class__.__name__)
        self.logger.setLevel(logging.DEBUG)

        logger_formatter = logging.Formatter("%(levelname)-10s:\t%(asctime)-20s:\t%(message)s", "%Y-%m-%d %H:%M:%S")
        if log_to_file:
            logger_file_handler = logging.FileHandler(os.path.join(self.log_dir, __class__.__name__ + '.log'), mode = 'w')
            logger_file_handler.setFormatter(logger_formatter)
            logger_file_handler.setLevel(logging.DEBUG)
            self.logger.addHandler(logger_file_handler)

        logger_stdout_handler = logging.StreamHandler(sys.stdout)
        logger_stdout_handler.setLevel(logging.INFO)

        self.logger.addHandler(logger_stdout_handler)

        self.df_uniprot_cols = ["Entry", "Uniprot_Sequence", "Query"]
        self.df_blast_col = ['blast_uniprot_id','blast_identity','fasta_id','mutated_Sequence','blast_seq','species','blast_fasta_id']

        self.blast_path = os.path.join(self.auxillary_dir,"ncbi-blast-2.13.0+","bin","blastp")
        self.uniprot_db = os.path.join(self.auxillary_dir,"UniprotKB-26042023")
        self._load_auxillary()


    def _load_auxillary(self):
        """
        load auxillary files: \'map_inchikey_to_canonicalSMILES.json\' and \'uniprot_sequences.csv\' and 
        put the result to attributes.

        Returns:
        --------
        map_inchikey_to_canonicalSMILES : dict
            loaded \'map_inchikey_to_canonicalSMILES.json\'.
        
        df_uniprot : pandas.DataFrame
            loaded \'uniprot_sequences.csv\'.
        """
        # map_inchikey_to_canonicalSMILES:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'r') as jsonfile:
                map_inchikey_to_canonicalSMILES = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_canonicalSMILES = {}

        # df_uniprot
        try:
            df_uniprot = pandas.read_csv(os.path.join(self.auxillary_dir, 'uniprot_sequences.csv'), sep = ';', index_col = None)
        except FileNotFoundError:
            df_uniprot = pandas.DataFrame([], columns = self.df_uniprot_cols)
            # raise FileNotFoundError('Support file: \'uniprot_sequences.csv\' not found. Please call self.update_df_uniprot or run uniprot_utils.py to download sequences.')

        assert len(df_uniprot.columns) == len(self.df_uniprot_cols)
        for i in range(len(df_uniprot.columns)):
            if df_uniprot.columns[i] != self.df_uniprot_cols[i]:
                raise ValueError('df_uniprot has different columns or column positions than in self.df_uniprot_cols')
        df_uniprot.set_index(self.df_uniprot_cols[0], drop = True, inplace = True)

        # df_blast
        try:
            df_blast = pandas.read_csv(os.path.join(self.auxillary_dir, 'df_blast.csv'), sep = ';', index_col = [0])
        except FileNotFoundError:
            df_blast = pandas.DataFrame([], columns = self.df_blast_col)

        assert len(df_blast.columns) == len(self.df_blast_col)
        for i in range(len(df_blast.columns)):
            if df_blast.columns[i] != self.df_blast_col[i]:
                raise ValueError('df_blast has different columns or column positions than in self.df_blast_cols')
        
        #
        self.map_inchikey_to_canonicalSMILES = map_inchikey_to_canonicalSMILES
        self.df_uniprot = df_uniprot
        self.df_blast = df_blast
        return map_inchikey_to_canonicalSMILES, df_uniprot, df_blast


    def _update_auxilary_map_inchikey_to_canonincalSMILES(self, candidate_idx):
        """
        update \'map_inchikey_to_canonicalSMILES.json\' with new InChI keys found in candidate_idx.

        Parameters:
        -----------
        candidate_idx : pandas.Index
            index with name \'InChI Key\' containing InChI keys to check for download. Info for the ones that are 
            missing in map_inchikey_to_canonicalSMILES.keys() are downloaded.
        """
        assert candidate_idx.name == 'InChI Key'
        map_inchikey_to_canonicalSMILES = self.map_inchikey_to_canonicalSMILES.copy()
        current_idx = pandas.Index(map_inchikey_to_canonicalSMILES.keys(), name = 'InChI Key')
        new_idx = candidate_idx.difference(current_idx)
        if len(new_idx) > 0:
            self.logger.info('Updating map_inchikey_to_canonicalSMILES...')
            NEW = get_map_inchikey_to_canonicalSMILES(new_idx.tolist())
            map_inchikey_to_canonicalSMILES.update(NEW)
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'w') as jsonfile:
                json.dump(map_inchikey_to_canonicalSMILES, jsonfile)
        return map_inchikey_to_canonicalSMILES


    def _update_auxillary_df_uniprot(self, full_df):
        """
        update \'uniprot_sequences.csv\' with new uniprot sequences found in full_df.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'Uniprot ID\' columns which will be used to check which IDs are new.
        
        Returns:
        --------
        df_uniprot : padnas.DataFrame
            updated df_uniprot
        """
        df_uniprot = self.df_uniprot.copy()
        candidate_idx = pandas.Index(full_df['Uniprot ID'].dropna().unique())
        new_idx = candidate_idx.difference(df_uniprot.index)
        if len(new_idx) > 0:
            self.logger.info('Updating df_uniprot...')
            NEW = get_uniprot_sequences(new_idx.tolist())
            NEW.set_index(self.df_uniprot_cols[0], drop = True, inplace = True)
            df_uniprot = df_uniprot.append(NEW, ignore_index = False, verify_integrity = True)
            df_uniprot.to_csv(os.path.join(self.auxillary_dir, 'uniprot_sequences.csv'), sep = ';', index = True)
        return df_uniprot

    def _update_auxillary_df_blast(self, full_df):
        """
        update \'df_blast.csv\' with new sequences found in full_df.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'mutated_Sequence\' columns which will be used to check which sequences are new.
        
        Returns:
        --------
        df_blast : padnas.DataFrame
            updated df_blast
        """
        df_blast = self.df_blast.copy()
        if df_blast.empty:
            self.logger.info('Creating df_blast...')
            new_sequences = full_df.drop_duplicates(subset=['species','mutated_Sequence']).copy()
            uniprot_db_path = [f"{p}/{s.replace(' ', '_')}/{s.replace(' ', '_')}" for p, s in zip([self.uniprot_db] * len(new_sequences.species.tolist()), new_sequences.species.tolist())]
            NEW = get_blast_data(new_sequences.mutated_Sequence.tolist(), uniprot_db_path, self.blast_path, new_sequences.species.tolist())
        else:
            all_sequences = full_df['mutated_Sequence'].dropna().unique().tolist()
            new_sequences = df_blast[~df_blast.mutated_Sequence.isin(all_sequences)]
            new_sequences.drop_duplicates(subset=['species','mutated_Sequence'], inplace=True)
            uniprot_db_path = [f"{p}/{s.replace(' ', '_')}/{s.replace(' ', '_')}" for p, s in zip([self.uniprot_db] * len(new_sequences.species.tolist()), new_sequences.species.tolist())]
            NEW = get_blast_data(new_sequences.mutated_Sequence.tolist(), uniprot_db_path, self.blast_path, new_sequences.species.tolist())
        if len(new_sequences) > 0:
            self.logger.info('Appending to df_blast...')
            df_blast = pandas.concat([df_blast, NEW], axis=0, ignore_index=True).reset_index(drop=True)
            df_blast.to_csv(os.path.join(self.auxillary_dir, 'df_blast.csv'), sep = ';', index = True)
        return df_blast

    def update_auxillary(self, df):
        """
        call and other operation necessery for all _update_auxillary_*.
        """
        _df = df.copy()
        _df = _df.dropna(subset = ['InChI Key'])
        candidate_idx = _df['InChI Key'].dropna()
        candidate_idx = candidate_idx.str.split(' ').explode() # TODO: pandas FutureWarning for this row.
        candidate_idx = pandas.Index(candidate_idx.unique(), name = 'InChI Key')
        self.map_inchikey_to_canonicalSMILES = self._update_auxilary_map_inchikey_to_canonincalSMILES(candidate_idx)
        self.df_uniprot = self._update_auxillary_df_uniprot(df)
        if 'mutated_Sequence' in df.columns:
            self.df_blast = self._update_auxillary_df_blast(df)
        return
    

    def add_implied_columns(self, df):
        """
        Add \'_Sequence\' column to df with retrieved sequences from UniProt (or keeping sequence if no UniProt ID is provided).

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        df_join : pandas.DataFrame
            dataframe with added column.
        """
        # self._load_auxillary()
        # self._update_auxillary(df)
        df_join = df.copy()
        if df_join["Uniprot ID"].isna().all():
            df_join['Uniprot_Sequence'] = float('nan')
            #print('empty')
        else:
            print(df_join.columns)
            df_join = df_join.join(self.df_uniprot[['Uniprot_Sequence']], on = 'Uniprot ID', how = 'left')
        
        df_join['_Sequence'] = df_join.apply(lambda x: merge_cols_with_priority(x, primary_col = 'Uniprot_Sequence', secondary_col = 'Sequence'), axis = 1)
        return df_join


    def get_map_inchikey_to_isomers(self, df):
        """
        For each InChI key that contains \'-UHFFFAOYSA-\' get all the isomers. 
        Isomers are identified using rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers

        If there are multiple InChI key in a record (i.e. mixture) then isomers are found for each InChI key separately.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe  to process

        Returns:
        --------
        map_inchikey_to_isomers : dict
            mapping from InChI key to list of isomers.
        """
        _df = df.copy()
        _df = _df.dropna(subset = ['InChI Key'])
        map_inchikey_to_isomers = {}
        for x in _df['InChI Key'].unique():
            _res = {}
            for inchikey in x.split(' '):
                if '-UHFFFAOYSA-' in inchikey: # This is done only for non-isomeric
                    _res[inchikey] = enumerate_isomers(self.map_inchikey_to_canonicalSMILES[inchikey])
            map_inchikey_to_isomers.update(_res)
        return map_inchikey_to_isomers


    def get_map_canonicalSMILES_to_isomers(self, df):
        """
        For each canonical SMILES get all the isomers. Isomers are identified 
        using rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers

        If there are multiple canonical SMILES in a record (i.e. mixture) then isomers are found for each SMILES separately.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe  to process

        Returns:
        --------
        map_canonicalSMILES_to_isomers : dict
            mapping from canonical SMILES to list of isomers.
        """
        _df = df.copy()
        _df = _df.dropna(subset = ['canonicalSMILES'])
        map_canonicalSMILES_to_isomers = {}
        for x in _df['canonicalSMILES'].unique():
            _res = {}
            for smiles in x.split(' '):
                _res[smiles] = enumerate_isomers(smiles)
            map_canonicalSMILES_to_isomers.update(_res)
        return map_canonicalSMILES_to_isomers


    @staticmethod
    def _update_mixture(x, _map_inchikey, _map_smiles):
        """
        Notes:
        ------
        NaN is kept and retured as NaN.
        """
        if x['InChI Key'] == x['InChI Key']:
            if ' ' in x['InChI Key']:
                return 'mixture'
            elif x['InChI Key'] in _map_inchikey.keys():
                if len(_map_inchikey[x['InChI Key']]) > 1:
                    return 'sum of isomers'
                else:
                    return 'mono'
            else:
                return 'mono'
        elif x['canonicalSMILES'] == x['canonicalSMILES']:
            if ' ' in x:
                return 'mixture'
            elif x['canonicalSMILES'] in _map_smiles.keys():
                if len(_map_smiles[x['canonicalSMILES']]) > 1:
                    return 'sum of isomers'
                else:
                    return 'mono'
            else:
                return 'mono'
        return float("nan")

    def update_mixture_col(self, df):
        """
        Create/update \'Mixture\' column based on isomers from map_inchikey_to_isomers and map_canonicalSMILES_to_isomers.

        If \'InChI Key\' (or \'canonicalSMILES\' if InChI key is not provided) contains space put label \'mixture\',
        else if number of isomers is greater than 1 put \'sum of isomers\' and 
        for other cases put \'mono\'.

        NaNs are kept as NaNs.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        new_df : pandas.DataFrame
            copy of df with updated \'Mixture\' column.
        """
        new_df = df.copy()
        map_inchikey_to_isomers = self.get_map_inchikey_to_isomers(df)
        map_canonicalSMILES_to_isomers = self.get_map_canonicalSMILES_to_isomers(df)
        new_df['Mixture'] = new_df.apply(lambda x: self._update_mixture(x, _map_inchikey = map_inchikey_to_isomers, _map_smiles = map_canonicalSMILES_to_isomers), axis = 1)
        return new_df   
        

    def update_mutated_sequence_col(self, df):
        """
        Create/update \'mutated_Sequence\' column based on Sequence, Uniprot ID and Mutation columns.

        First _Sequence column is populated with sequences from Uniprot and from the data and then the mutation 
        is performed on these sequences.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        new_df : pandas.DataFrame
            copy of df with updated \'Mixture\' column.
        """
        new_df = df.copy()
        new_df_with_implied = self.add_implied_columns(df)
        new_df['mutated_Sequence'] = new_df_with_implied.apply(lambda x: perform_mutation(x, mutation_col = 'Mutation', seq_col = '_Sequence'), axis = 1)
        return new_df

    def add_blast_data(self, df):
        """
        Create/update \'mutated_Sequence\' column based on Sequence, Uniprot ID and Mutation columns.

        First _Sequence column is populated with sequences from Uniprot and from the data and then the mutation 
        is performed on these sequences.

        Paramters:
        ----------
        df : pandas.DataFrame
            dataframe being processed.

        Returns:
        --------
        new_df : pandas.DataFrame
            copy of df with updated \'Mixture\' column.
        """
        new_df = df.copy()
        before = len(new_df)
        new_df = new_df.merge(self.df_blast, on=["mutated_Sequence","species"], how="inner")
        if before != len(new_df):
            print(before, len(new_df))
            raise ValueError 
        return new_df


    def __call__(self, df):
        """
        Paramters:
        ----------
        df : pandas.DataFrame
            raw dataframe.
    
        Returns:
        --------
        new_df : pandas.DataFrame
            formatted dataframe
        """
        self.update_auxillary(df)
        new_df = self.update_mixture_col(df)
        new_df = self.update_mutated_sequence_col(new_df)
        self.update_auxillary(new_df)
        new_df = self.add_blast_data(new_df)
        return new_df
            





if __name__ == '__main__':
    df = pandas.read_csv('./RawData/bc_Max/Clean/Araneda_2004_done_20220221.csv', sep=';', index_col = 0)

    formatter = PreFormatter()
    test, exclude_test = formatter(df)

    post_formatter = PostFormatter()
    new_df = post_formatter(test)
    print(new_df)