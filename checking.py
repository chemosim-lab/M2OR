import os
import sys
import numpy
import pandas
import logging
import re
import json

from utils import perform_mutation, merge_cols_with_priority, enumerate_isomers
from uniprot_utils import get_uniprot_sequences
from pubchem_utils import get_map_inchikey_to_CID, get_map_inchikey_to_synonyms, get_map_name_to_inchikeys, get_map_inchikey_to_canonicalSMILES

_logging_file_path = 'Log file path'



# Create logger


# (OK) TODO: Check consistency : group by mutated_seq and check len(unique(mutated_ID)) , group by mutated_ID and check len(unique(mutated_seq))
# (OK) TODO: Mixture: mono + multiple keys
# (OK) TODO: Same everything except Response.
# (OK) TODO: Inchi Key not on pubchem.
# (OK) TODO: check pubchempy for going from inchikey to name and back. It can happen that we have incorrect inchi key ! ! ! !


class Checker:
    """
    class running all the consistency checking before merging with clean data. The main idea of checker is to highlight discrepancies and then users
    needs to fix them themselves.

    Each check is implemented as a separate method and can in principle be used separately (although there are specific asumptions 
    for each check like non-NaN values in entries.)

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

    map_inchikey_to_CID_cols : list
        columns expected to be found in \'map_inchikey_to_CID.csv\'. This serves as a precaution.

    logging_cols : list
        list of columns to include in examples in logging output.

    not_nan_cols : list
        list of column names that should not contain NaN.

    conditioned_not_nan_cols : list
        list of dictionaries with structure {'cond_col' : colname, 'cond_val' : value, 'cols' : [list_of_column_names]} that
        defines condition: \"For all records that have \'cond_val\' in \'cond_col\' check if entries in \'cols\' are not NaN.\"

    format_cols : list
        list of dictionaries with structure {'col' : colname, 'pattern' : regex_pattern, 'sep' : multiple_entries_separator} that
        are used to check the formatting in \'col\'. If there are multiple elements in one entry (e.g. multiple InChI keys for mixtures)
        separator \'sep\' is used to split entry into elements.

    castable_cols : list
        list of dictionaries with structure {'col' : colname, 'Type' : data_type, 'except_values' : values_to_ignore, 'ignore_patterns' : patterns_to_replace_with_EMPTY, 'sep' : separator_for_split}
        that are used to check if entries in a given column can be casted to a data type \'Type\'. Exception of the rule are given by \'except_values\' and special patterns like \'>\' can be ignored 
        using \'ignore_patterns\'. Use None for \'except_values\' and \'ignore_patterns\' if there is no exception.

    not_castable_cols : list
        list of dictionaries with structure {'col' : colname, 'Type' : data_type, 'except_values' : values_to_ignore, 'ignore_patterns' : patterns_to_replace_with_EMPTY, 'sep' : separator_for_split} 
        that is used to check if entries in a given column are not of a given data type \'Type\'. Exceptions of the rule are given by \'except_values\' and special patterns like can be ignored 
        using \'ignore_patterns\'. Use None for \'except_values\' and \'ignore_patterns\' if there is no exception.

    df_uniprot : pandas.DataFrame
        auxillary dataframe with mapping from uniprot ID to sequence. It corresponds to \'uniprot_sequences.csv\'.

    map_inchikey_to_CID : pandas.DataFrame
        auxillary dataframe with mapping from InChI key to CID. It corresponds to \'map_inchikey_to_CID.csv\'.

    map_inchikey_to_canonicalSMILES : dict
        auxillary dictionary mapping InChI key to canonical SMILES.

    map_inchikey_to_synonyms : dict
        auxillary dictionary mapping InChI key to synonyms for pubchem.

    map_name_to_inchikeys : dict
        auxillary dictionary mapping name to InChI key. NOTE: Not used now

    check_results : dict
        dictionary with the results of the performed checks.

    References:
    -----------
    InChI Key format: https://gist.github.com/lsauer/1312860/264ae813c2bd2c27a769d261c8c6b38da34e22fb
    """
    def __init__(self, auxillary_dir = 'Data', log_dir = 'logs'):
        self.auxillary_dir = auxillary_dir
        self.log_dir = log_dir

        self._init_logger(__class__.__name__)
        self._init_config()
        self.check_results = {}


    def _init_logger(self, logger_name):
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(logging.DEBUG)

        logger_formatter = logging.Formatter("%(levelname)-10s:\t%(asctime)-20s:\t%(message)s", "%Y-%m-%d %H:%M:%S")
        logger_file_handler = logging.FileHandler(os.path.join(self.log_dir, logger_name + '.log'), mode = 'w')
        logger_file_handler.setFormatter(logger_formatter)
        logger_stdout_handler = logging.StreamHandler(sys.stdout)

        logger_file_handler.setLevel(logging.DEBUG)
        logger_stdout_handler.setLevel(logging.INFO)

        self.logger.addHandler(logger_file_handler)
        self.logger.addHandler(logger_stdout_handler)


    def _init_config(self):
        """
        create static attributes.
        """
        self.df_uniprot_cols = ["Entry", "Uniprot_Sequence", "Query"]
        self.map_inchikey_to_CID_cols = ['InChI Key', 'CID']
        # self.map_name_to_inchikey_cols = ['Name', 'InChI Key']

        self.logging_cols = ['species',
                            'Mutation',
                            'Gene ID',
                            'Uniprot ID',
                            # 'Sequence',
                            'Name',
                            # 'CID',
                            # 'CAS',
                            'InChI Key',
                            'Parameter',
                            'Value',
                            'Unit',
                            'Value_Screen',
                            'Unit_Screen',
                            'Responsive',
                            # 'Type',
                            'Reference',
                            'DOI',
                            'Reference Position',
                            'Assay',
                            'Delivery',
                            'Tag',
                            'Gprotein',
                            'Cell_line'
                            #'Mixture',
                            ]
        
        # Config:
        self.not_nan_cols = ['_MolID', 
                            'Parameter', 
                            'Responsive', 
                            'DOI',
 #                           'Mixture',
                            '_Sequence',
                            'nbr_measurements',
                            'Type',
                            'Assay',
                            'Delivery',
                            'Cell_line'
                            ]

        # For a given conditioning value cond_val in a given 
        # conditioning column cond_col the values in cols must NOT be NaN:
        self.conditioned_not_nan_cols = [{'cond_col' : 'Parameter', 'cond_val' : 'ec50',    'cols' : ['Value', 'Unit']},
                                    {'cond_col' : 'Parameter', 'cond_val' : 'screening',    'cols' : ['Value_Screen', 'Unit_Screen']},
                                    {'cond_col' : 'Parameter', 'cond_val' : 'rlu/rlumax',   'cols' : ['Value_Screen', 'Unit_Screen']},
                                    {'cond_col' : 'Parameter', 'cond_val' : 'potential',    'cols' : ['Value', 'Unit', 'Value_Screen', 'Unit_Screen']},
                                    {'cond_col' : 'Parameter', 'cond_val' : 'calcium signal',   'cols' : ['Value_Screen', 'Unit_Screen']},
                                    ]

        self.format_cols = [{'col' : 'InChI Key', 'pattern' : r'^([A-Z\-]{27})$', 'sep' : ' '},
                            {'col' : 'Mutation',  'pattern' : r'^([AGILPVFWYDERHKSTCMNQ][0-9]{1,3}[AGILPVFWYDERHKSTCMNQ])$', 'sep' : '_'},
                            {'col' : 'Sequence',  'pattern' : r'[AGILPVFWYDERHKSTCMNQ]+', 'sep' : ' '},
                        ]

        self.castable_cols = [{'col' : 'Value',         'Type' : float, 'except_values' : ['n.d'],  'ignore_patterns' : ['>', '<', '-', '.'],   'sep' : ' '}, # TODO: Check ignore_patterns '.' ..it is here because of '0.1-0.3'. Is there a better way to treat this?
                              {'col' : 'Value_Screen',  'Type' : float, 'except_values' : None,     'ignore_patterns' : None,                   'sep' : ' '},
                              {'col' : 'Responsive',    'Type' : int,   'except_values' : None,     'ignore_patterns' : None,                   'sep' : None},
                              {'col' : 'CID',           'Type' : int,   'except_values' : None,     'ignore_patterns' : None,                   'sep' : None},
                        ]

        self.not_castable_cols = [{'col' : 'Parameter',    'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'Unit',         'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'Unit_Screen',  'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'Type',         'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'Reference',    'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'DOI',          'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                  {'col' : 'Mixture',      'Type' : float, 'except_values' : None, 'ignore_patterns' : None, 'sep' : None},
                                ]

        # String for some columns should be within predefined categories
        self.categorical_values = [{'col' : 'Type',      'allowed_values':['Ca2+','Luc','cAMP','SEAP','I','Conductance','GFP']},
                                   {'col' : 'Cell_line', 'allowed_values':['HEK','H3A','Ocy','OSN','OB','HeLa/Olf','NxG108CC15','ScL21','Yeast','HEPG2','HUH7','LNCAP']},
                                   {'col' : 'Assay',     'allowed_values':['in vitro','in vivo','ex vivo']},
                                   {'col' : 'Delivery',  'allowed_values':['liquid','gaz']},
                                   {'col' : 'Gprotein',  'allowed_values':['olf','gα15/gα16','gα16','gαq']},
                                   {'col' : 'Tag',       'allowed_values':['rho','myc','flag','rho lucy','gfp','il-6-halotag']},
                                   {'col' : 'Parameter', 'allowed_values':['ec50','raw','norm_other','norm_pair','norm_rec','norm_mol']}
                                ]


    def _logging_format_dataframe(self, df):
        """
        transform fail examples dataframe to text for logging.
        """
        n_fails = len(df)
        msg = 'num. of fails: {}, fail examples: \n' + df[self.logging_cols].head(20).to_string(max_colwidth = 50)
        return msg.format(n_fails)


    def _load_auxillary(self):
        """
        load auxillary files: \'uniprot_sequences.csv\', \'map_inchikey_to_CID.csv\', \'map_inchikey_to_canonicalSMILES.json\', \'map_inchikey_to_synonyms.json\' 
        and put the result to attributes.

        Returns:
        --------
        df_uniprot : pandas.DataFrame
            loaded  \'uniprot_sequences.csv\'
        
        map_inchikey_to_CID : pandas.DataFrame
            loaded  \'map_inchikey_to_CID.csv\'
        
        map_inchikey_to_canonicalSMILES : dict
            loaded \'map_inchikey_to_canonicalSMILES.json\'.

        map_inchikey_to_synonyms : dict
            loaded \'map_inchikey_to_synonyms.json\'
        """
        # df_uniprot:
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

        # map_inchikey_to_CID:
        try:
            map_inchikey_to_CID = pandas.read_csv(os.path.join(self.auxillary_dir, 'map_inchikey_to_CID.csv'), sep = ';', index_col = None)
        except FileNotFoundError:
            map_inchikey_to_CID = pandas.DataFrame([], columns = self.map_inchikey_to_CID_cols)
            # raise FileNotFoundError('Support file: \'map_inchikey_to_CID.csv\' not found. Please run pubchem_utils.py to download inchi key to cid mapping.')

        assert len(map_inchikey_to_CID.columns) == len(self.map_inchikey_to_CID_cols)
        for i in range(len(map_inchikey_to_CID.columns)):
            if map_inchikey_to_CID.columns[i] != self.map_inchikey_to_CID_cols[i]:
                raise ValueError('map_inchikey_to_CID has different columns or column positions than in self.map_inchikey_to_CID_cols')
        map_inchikey_to_CID.set_index(self.map_inchikey_to_CID_cols[0], drop = True, inplace = True)
        map_inchikey_to_CID = map_inchikey_to_CID.squeeze()

        # map_inchikey_to_canonicalSMILES:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'r') as jsonfile:
                map_inchikey_to_canonicalSMILES = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_canonicalSMILES = {}

        # map_inchikey_to_synonyms:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_synonyms.json'), 'r') as jsonfile:
                map_inchikey_to_synonyms = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_synonyms = {}

        # map_name_to_inchikeys:
        # try:
        #     with open(os.path.join(self.auxillary_dir, 'map_name_to_inchikeys.json'), 'r') as jsonfile:
        #         map_name_to_inchikeys = json.load(jsonfile)
        # except FileNotFoundError:
        #     map_name_to_inchikeys = {}

        #
        self.df_uniprot = df_uniprot
        self.map_inchikey_to_CID = map_inchikey_to_CID
        self.map_inchikey_to_canonicalSMILES = map_inchikey_to_canonicalSMILES
        self.map_inchikey_to_synonyms = map_inchikey_to_synonyms
        # self.map_name_to_inchikeys = map_name_to_inchikeys
        return df_uniprot, map_inchikey_to_CID, map_inchikey_to_canonicalSMILES, map_inchikey_to_synonyms # , map_name_to_inchikeys

    def _update_auxilary_df_uniprot(self, full_df):
        """
        update \'uniprot_sequences.csv\' with new uniprot sequences found in full_df.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'Uniprot ID\' column which will be used to check which IDs are new.

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

    def _update_auxilary_map_inchikey_to_CID(self, full_df):
        """
        update \'map_inchikey_to_CID.csv\' with new InChI Keys found in full_df.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'InChI Key\' column which will be used to check which IDs are new.

        Returns:
        --------
        map_inchikey_to_CID : padnas.DataFrame
            updated map_inchikey_to_CID
        """
        map_inchikey_to_CID = self.map_inchikey_to_CID.copy()
        candidate_idx = full_df['InChI Key'].dropna()
        candidate_idx = candidate_idx.str.split(' ').explode() # TODO: pandas FutureWarning for this row.
        candidate_idx = pandas.Index(candidate_idx.unique())
        new_idx = candidate_idx.difference(map_inchikey_to_CID.index)
        if len(new_idx) > 0:
            self.logger.info('Updating map_inchikey_to_CID...')
            NEW = get_map_inchikey_to_CID(new_idx.tolist())
            NEW = pandas.Series(NEW, dtype = float)
            NEW.index.name = self.map_inchikey_to_CID_cols[0] # InChI Key
            NEW.name = self.map_inchikey_to_CID_cols[1] # CID
            map_inchikey_to_CID = map_inchikey_to_CID.append(NEW, ignore_index = False, verify_integrity = True)
            map_inchikey_to_CID.to_csv(os.path.join(self.auxillary_dir, 'map_inchikey_to_CID.csv'), sep = ';')
        return map_inchikey_to_CID

    def _update_auxilary_map_inchikey_to_canonicalSMILES(self, full_df):
        """
        update \'map_inchikey_to_canonicalSMILES.json\' with new InChI keys found in candidate_idx.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'InChI Key\' column which will be used to check which IDs are new.
        
        Returns:
        --------
        map_inchikey_to_canonicalSMILES : dict
            updated map_inchikey_to_canonicalSMILES
        """
        map_inchikey_to_canonicalSMILES = self.map_inchikey_to_canonicalSMILES.copy()
        current_idx = pandas.Index(map_inchikey_to_canonicalSMILES.keys(), name = 'InChI Key')
        candidate_idx = full_df['InChI Key'].dropna()
        candidate_idx = candidate_idx.str.split(' ').explode() # TODO: pandas FutureWarning for this row.
        candidate_idx = pandas.Index(candidate_idx.unique())
        new_idx = candidate_idx.difference(current_idx)
        if len(new_idx) > 0:
            self.logger.info('Updating map_inchikey_to_canonicalSMILES...')
            NEW = get_map_inchikey_to_canonicalSMILES(new_idx.tolist())
            map_inchikey_to_canonicalSMILES.update(NEW)
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'w') as jsonfile:
                json.dump(map_inchikey_to_canonicalSMILES, jsonfile)
        return map_inchikey_to_canonicalSMILES

    def _update_auxilary_map_inchikey_to_synonyms(self, full_df):
        """
        update \'map_inchikey_to_synonyms.json\' with new InChI keys found in candidate_idx.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe that is being processed. It needs to have \'InChI Key\' column which will be used to check which IDs are new.
        
        Returns:
        --------
        map_inchikey_to_synonyms : dict
            updated map_inchikey_to_synonyms
        """
        map_inchikey_to_synonyms = self.map_inchikey_to_synonyms.copy()
        current_idx = pandas.Index(map_inchikey_to_synonyms.keys(), name = 'InChI Key')
        candidate_idx = full_df['InChI Key'].dropna()
        candidate_idx = candidate_idx.str.split(' ').explode() # TODO: pandas FutureWarning for this row.
        candidate_idx = pandas.Index(candidate_idx.unique())
        new_idx = candidate_idx.difference(current_idx)
        if len(new_idx) > 0:
            self.logger.info('Updating map_inchikey_to_synonyms...')
            NEW = get_map_inchikey_to_synonyms(new_idx.tolist())
            map_inchikey_to_synonyms.update(NEW)
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_synonyms.json'), 'w') as jsonfile:
                json.dump(map_inchikey_to_synonyms, jsonfile)
        return map_inchikey_to_synonyms

    # def _update_auxilary_map_name_to_inchikeys(self, full_df):
    #     map_name_to_inchikeys = self.map_name_to_inchikeys.copy()
    #     current_idx = pandas.Index(map_name_to_inchikeys.keys(), name = 'Name')
    #     candidate_idx = full_df[full_df['Mixture'] == 'mono']['Name'].dropna() # TODO: Can we somehow work with mixtures?
    #     candidate_idx = pandas.Index(candidate_idx.unique())
    #     new_idx = candidate_idx.difference(current_idx)
    #     if len(new_idx) > 0:
    #         self.logger.info('Updating map_name_to_inchikeys...')
    #         NEW = get_map_name_to_inchikeys(new_idx.tolist())
    #         map_name_to_inchikeys.update(NEW)
    #         with open(os.path.join(self.auxillary_dir, 'map_name_to_inchikeys.json'), 'w') as jsonfile:
    #             json.dump(map_name_to_inchikeys, jsonfile)
    #     return map_name_to_inchikeys

    
    def _update_auxillary(self, full_df):
        """
        call and other operation necessery for all _update_auxillary_*.
        """
        self.df_uniprot = self._update_auxilary_df_uniprot(full_df)
        self.map_inchikey_to_CID = self._update_auxilary_map_inchikey_to_CID(full_df)
        self.map_inchikey_to_canonicalSMILES = self._update_auxilary_map_inchikey_to_canonicalSMILES(full_df)
        self.map_inchikey_to_synonyms = self._update_auxilary_map_inchikey_to_synonyms(full_df)
        # self.map_name_to_inchikey = self._update_auxilary_map_name_to_inchikeys(full_df)
        return 


    def add_implied_columns(self, full_df):
        """
        Add \'_Sequence\' and \'_MolID\' columns to df. 
        
        \'_Sequence\' column has retrieved sequences from UniProt in it (or keeping sequence if no UniProt ID is provided). 
        \'_MolID\' column has InChI Key if available and canonical SMILES if InChI key is missing.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        df_join : pandas.DataFrame
            dataframe with added columns.
        """
        self._load_auxillary()
        self._update_auxillary(full_df)
        df_join = full_df.copy()
        if df_join["Uniprot ID"].isna().all():
            df_join['Uniprot_Sequence'] = float('nan')
        else:
            df_join = df_join.join(self.df_uniprot[['Uniprot_Sequence']], on = 'Uniprot ID', how = 'left')
        df_join['_Sequence'] = df_join.apply(lambda x: merge_cols_with_priority(x, primary_col = 'Uniprot_Sequence', secondary_col = 'Sequence'), axis = 1)
        df_join['_MolID'] = df_join.apply(lambda x: merge_cols_with_priority(x, primary_col = 'InChI Key', secondary_col = 'canonicalSMILES'), axis = 1)
        return df_join


    def check_not_nan(self, full_df):
        """
        Check if there are any missing values in columns that should always be filled.

        Non-NaN columns are defined in \'not_nan_cols\'.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_not_nan')
        passed = True
        clean_df = full_df.copy()
        for case in self.not_nan_cols:
            _df = clean_df
            condition = _df[case].isna()
            if condition.any():
                passed = False
                fail_example = _df[condition]
                self.logger.debug('FAIL in check_not_nan:  col: \'{}\''.format(case))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)] # Clean from found errors.
        if passed:
            self.logger.info('FINISHED: check_not_nan:  PASS')
        else:
            self.logger.warning('FINISHED: check_not_nan:  FAIL')
        return clean_df, passed


    def check_conditioned_not_nan(self, full_df):
        """
        For each entry in \'conditioned_not_nan_cols\', for records where value in \'cond_col\' is \'cond_val\' check if value in other 
        columns \'cols\' is missing.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_conditioned_not_nan')
        passed = True
        clean_df = full_df.copy()
        for case in self.conditioned_not_nan_cols:
            _df = clean_df[clean_df[case['cond_col']] == case['cond_val']]
            condition = _df[case['cols']].isna().any(axis = 1)
            if condition.any():
                passed = False
                fail_example = _df[condition]
                self.logger.debug('FAIL in check_conditioned_not_nan:  cond_col: \'{}\'  cond_val: \'{}\'  cols:  \'{}\''.format(*case.values()))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_conditioned_not_nan:  PASS')
        else:
            self.logger.warning('FINISHED: check_conditioned_not_nan:  FAIL')
        return clean_df, passed


    @staticmethod
    def _check_castable(x, _type):
        """
        Try to cast x to _type.
        """
        try:
            _type(x)
            return True
        except ValueError:
            return False

    def _check_castable_multiple(self, x, _type, sep = None):
        """
        Try if input x is castable to a given data type _type. If sep is not None, x is first split and then each element is checked separately.
        x passes the check if all elements inside are passing _check_castable.

        This function is expected to be used in pandas apply.

        Parameters:
        -----------
        x : Any
            input

        _type : object
            data type to try to cast to.

        sep : str, optional (default=None)
            separator if there are multiple entries.

        Returns:
        --------
        passed : bool
            result of the check.
        """
        passed = True
        if isinstance(x, str) and sep is not None:
            for ele in x.split(sep):
                if not self._check_castable(ele, _type):
                    passed = False
        else:
            if not self._check_castable(x, _type):
                passed = False
        return passed

    def check_castable(self, full_df):
        """
        Check if entries in some columns are castable to a given data type. NaNs are ignored.

        For each column in \'castable_cols\', first take out exceptions from the rule (given by \'except_values\') and replace patterns to ignore (\'ignore_patterns\') with 
        empty string. Then check if entries in the column are castable to a given data type \'Type\'. If there are multiple values (e.g. for mixtures) castable check is
        performed for all of them separately. They are split using \'sep\'.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_castable')
        passed = True
        clean_df = full_df
        for case in self.castable_cols:
            _df = clean_df
            _df = _df.dropna(subset = [case['col']])
            if case['except_values'] is not None:
                for except_value in case['except_values']:
                    _df = _df[~(_df[case['col']] == except_value)]
            if case['ignore_patterns'] is not None:
                for pat in case['ignore_patterns']:
                    _df[case['col']] = _df[case['col']].astype(str).str.replace(pat, '', regex = True)
            # if case['col'] == "Value_Screen":
            #     condition = _df[case['col']].apply(lambda x: x.replace(" ","").replace('.','').isdigit() if isinstance(x, str) else isinstance(x, (int, float))) #Split and check (see mutation)
            # else:
            #     condition = _df[case['col']].apply(lambda x: self._check_castable(x, _type = case['Type']))
            condition = _df[case['col']].apply(lambda x: self._check_castable_multiple(x, _type = case['Type'], sep = case['sep']))
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_castable:  col: \'{}\'  type: \'{}\''.format(*case.values()))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_castable:  PASS')
        else:
            self.logger.warning('FINISHED: check_castable:  FAIL')
        return clean_df, passed


    def _check_not_castable_multiple(self, x, _type, sep = None):
        """
        Try if input x is castable to a given data type _type and if yes the check fails. If sep is not None, x is first split and then each 
        element is checked separately. x passes the check if all elements inside are *not* passing _check_castable.

        This function is expected to be used in pandas apply.

        Parameters:
        -----------
        x : Any
            input

        _type : object
            data type to try to cast to.

        sep : str, optional (default=None)
            separator if there are multiple entries.

        Returns:
        --------
        passed : bool
            result of the check.
        """
        passed = True
        if isinstance(x, str) and sep is not None:
            for ele in x.split(sep):
                if self._check_castable(ele, _type):
                    passed = False
        else:
            if self._check_castable(x, _type):
                passed = False
        return passed


    def check_not_castable(self, full_df):
        """
        Check if entries in some columns are not of a given data type. NaNs are ignored.

        For each column in \'not_castable_cols\', first take out exceptions from the rule (given by \'except_values\') and replace patterns to ignore (\'ignore_patterns\') with 
        empty string. Then check if entries in the column are castable to a given data type \'Type\' and if yes then check is *not* passed. If there are multiple values 
        (e.g. for mixtures) this check is performed for all of them separately. They are split using \'sep\'.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_not_castable')
        passed = True
        clean_df = full_df
        for case in self.not_castable_cols:
            _df = clean_df
            _df = _df.dropna(subset = [case['col']]) # TODO: Is this necessary since we clean df along the way (NaN can be tested before this)? ?
            if case['except_values'] is not None:
                for except_value in case['except_values']:
                    _df = _df[~(_df[case['col']] == except_value)]
            if case['ignore_patterns'] is not None:
                for pat in case['ignore_patterns']:
                    _df[case['col']] = _df[case['col']].str.replace(pat, '', regex = True)
            # condition = ~_df[case['col']].apply(lambda x: self._check_castable(x, _type = case['Type']))
            condition = _df[case['col']].apply(lambda x: self._check_not_castable_multiple(x, _type = case['Type'], sep = case['sep']))
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_not_castable:  col: \'{}\'  type: \'{}\''.format(*case.values()))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_not_castable:  PASS')
        else:
            self.logger.warning('FINISHED: check_not_castable:  FAIL')
        return clean_df, passed

    
    @staticmethod
    def _check_concat_format(x, pat, sep):
        if sep is None:
            sep = ' '
        return all([True if re.match(pat, s) else False for s in x.split(sep)])

    def check_format(self, full_df):
        """
        Check if entries in columns given by \'format_cols\' have correct format. For example if InChI Key is 27 characters 
        long with only letters and three dashes (see \'format_cols\' for precise format). If there are multiple elements in one entry,
        each element is check separately. Elements are separated by \'sep\'. This function uses regular expression.
        NaNs are ignored.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_format')
        passed = True
        clean_df = full_df
        for case in self.format_cols:
            _df = clean_df
            _df = _df.dropna(subset = [case['col']]) # NaNs are checked in check_not_nan
            condition = _df[case['col']].apply(lambda x: self._check_concat_format(x, case['pattern'], case['sep']))
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_format:  col: \'{}\'  pattern: \'{}\'  sep:  \'{}\''.format(*case.values()))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_format:  PASS')
        else:
            self.logger.warning('FINISHED: check_format:  FAIL')
        return clean_df, passed


    def check_mixture_format(self, full_df):
        """
        Check if there is no separator in records that are \'mono\' or \'sum of isomers\'.
        Check if there is a separator in records that are \'mixture\'.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_mixture_format')
        passed = True
        clean_df = full_df

        sep = ' ' # TODO check separator here

        _df = clean_df[clean_df['Mixture'].str.lower() != 'mixture']
        # condition_mono_sum = ~_df['InChI Key'].dropna().str.contains(sep)
        condition_mono_sum = ~_df['_MolID'].dropna().str.contains(sep)
        if not condition_mono_sum.all():
            passed = False
            self.logger.debug('FAIL in check_mixture_format: Mixture: \'mono or sum of isomers\', \'_MolID\' contains \'{}\''.format(sep))
            fail_example = _df[~condition_mono_sum]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        _df = clean_df[clean_df['Mixture'].str.lower() == 'mixture']
        # condition_mix = _df['InChI Key'].dropna().str.contains(sep)
        condition_mix = _df['_MolID'].dropna().str.contains(sep)
        if not condition_mix.all():
            passed = False
            self.logger.debug('FAIL in check_mixture_format: Mixture:  \'mixture\', \'InChI Key\' without \'{}\''.format(sep))
            fail_example = _df[~condition_mix]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_mixture_format:  PASS')
        else:
            self.logger.warning('FINISHED: check_mixture_format:  FAIL')
        return clean_df, passed

    @staticmethod
    def _check_inchikey_on_pubchem(x, list_inchikeys):
        passed = True
        if x == x:
            x = x.split()
            for inchikey in x:
                if inchikey not in list_inchikeys:
                    passed = False
        return passed

    def check_inchikey_on_pubchem(self, full_df):
        """
        Check if InChI keys in \'InChI Key\' column can be found on PubChem. This is using \'map_inchikey_to_CID\' and 
        internally pubchempy.get_compounds (see get_map_inchikey_to_CID in pubchem_utils.py). For multiple elements the
        check is performed on each separately. Separator is assumed to be *space*.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_inchikey_on_pubchem')
        passed = True
        clean_df = full_df.copy()
        _df = clean_df
        condition = _df['InChI Key'].apply(lambda x: self._check_inchikey_on_pubchem(x, list(self.map_inchikey_to_CID.index)))
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_inchikey_on_pubchem:  PASS')
        else:
            self.logger.warning('FINISHED: check_inchikey_on_pubchem:  FAIL')
        return clean_df, passed


    @staticmethod
    def _check_chirality(x, _map):
        """
        Check if molecules that are supposed to be achiral have more than 1 stereoisomer. 

        InChI keys containing \'-UHFFFAOYSA-\' should be achiral as well as all having cacnonical SMILES.

        Parameters:
        -----------
        x : pandas.Series
            row of the dataframe.

        _map : dict
            mapping from InChI key to canonical SMILES.
        """
        passed = True
        if x['InChI Key'] == x['InChI Key']:
            for inchikey in x['InChI Key'].split(' '):
                if '-UHFFFAOYSA-' in inchikey:
                    isomers = enumerate_isomers(_map[inchikey])
                    # mol = Chem.MolFromSmiles(_map[inchikey])
                    # chiralCenters = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True, includeCIP=False, useLegacyImplementation=False)
                    # isomers = tuple(EnumerateStereoisomers(mol))
                    if len(isomers) > 1:
                        passed = False
        elif x["canonicalSMILES"] == x["canonicalSMILES"]:
            isomers = enumerate_isomers(x["canonicalSMILES"])
            if len(isomers) > 1:
                passed = False
        return passed

    def check_chirality(self, full_df):
        """
        For InChI keys containing \'-UHFFFAOYSA-\' and for canonical SMILES check if number of stereoisomers is 1. These records should not have isomers
        and check fails if there are any (i.e. more than 1). For multiple elements the check is performed for each element and separator is assumed to
        be *space*. NaNs are ignored. Isomers are identified using rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers (see enumerate_isomers is utils.py).

        InChI keys are mapped to canonical SMILES using \'map_inchikey_to_canonicalSMILES\' and stereoisomers are identified using these SMILES.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.

        Notes:
        ------
        Doesn't work on mixtures.
        """
        self.logger.info('STARTED: check_chirality')
        passed = True
        clean_df = full_df.copy()

        _df = clean_df[clean_df['Mixture'] == "mono"]
        condition = _df.apply(lambda x: self._check_chirality(x, _map = self.map_inchikey_to_canonicalSMILES), axis=1)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        _df = clean_df[clean_df['Mixture'] == "sum of isomers"]
        condition = ~_df.apply(lambda x: self._check_chirality(x, _map = self.map_inchikey_to_canonicalSMILES), axis=1)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        
        if passed:
            self.logger.info('FINISHED: check_chirality:  PASS')
        else:
            self.logger.warning('FINISHED: check_chirality:  FAIL')
        return clean_df, passed

    @staticmethod
    def _clean_string_name(x):
        """
        Perform string cleaning and normalization on molecule's name.
        """
        try:
            x_clean = re.sub(r'\s+','', x.strip()).lower()
            x_clean = re.sub(r'[0-9]{0,1}r\/[0-9]{0,1}s', '', x_clean)
            x_clean = re.sub(r'[0-9]{0,1}e,[0-9]{0,1}z', '', x_clean)
            x_clean = re.sub(r'[0-9]{1}[ez]', '', x_clean)
            x_clean = x_clean.replace(",sumofisomers", "").replace("(+/-)-","").replace("-","").replace("+","")
            x_clean = x_clean.replace("\u03b1","alpha").replace("\u03B4","delta").replace("\u03B3","gamma").replace("\u03B2","beta")
            x_clean = x_clean.replace("d","").replace("l","")
            x_clean = x_clean.replace("(","").replace(")","")
        except Exception as e:
            print(x)
            raise e
        return x_clean

    def _check_inchikey_vs_name(self, x, cond_col, test_col, _map):
        """
        Check if text in \'test_col\' is in the list of possible synonyms for \'cond_col\'. It is mainly used to check if name of molecule is in synonyms retrieved
        from PubChem by InChI key.

        This is expected to be used with pandas apply.

        Parameters:
        -----------
        x : pandas.Series
            row of a dataframe

        cond_col : str
            name of the column containing identifier. Entry from here are mapped to synonyms using _map.

        test_col : str
            name of the column with text that should be in the list of synonyms.

        _map : dict
            mapping form identifier to synonyms.

        Notes:
        ------
        This works only with non-mixtures only, because mixtures can have name that is unrelated to names of the elements inside.
        *Entries in \'cond_col\' that are not in _map.keys() are ignored*.
        """
        # TODO: This is not checking for cols not in map.keys() ! ! !
        passed = True
        if x[cond_col] in _map.keys():
            x_clean = self._clean_string_name(x[test_col])
            
            for synonym in _map[x[cond_col]]:
                try: 
                    self._clean_string_name(synonym)
                except Exception as e:
                    print(_map[x[cond_col]], x[cond_col])
                    raise e
            if not x_clean in [self._clean_string_name(synonym) for synonym in _map[x[cond_col]]]:
                passed = False
        return passed

    def check_inchikey_vs_name(self, full_df):
        """
        Check if names can be found inside synonyms retrieved by InChI Key.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.

        Notes:
        ------
        This works only with non-mixtures only, because mixtures can have name that is unrelated to names of the elements inside.
        """
        # TODO: For now this is checked only for non-mixtures.
        self.logger.info('STARTED: check_inchikey_vs_name')
        passed = True
        clean_df = full_df.copy()

        # _df = clean_df[clean_df['Mixture'] == 'mixture']
        # _df = _df.dropna(subset = ['Name'])
        # condition_name = _df.apply(lambda x: self._check_inchikey_vs_name(x, cond_col = 'Name', test_col = 'InChI Key', _map = self.map_name_to_inchikeys), axis = 1)
        # if not condition_name.all():
        #     passed = False
        #     self.logger.debug('FAIL in check_inchikey_vs_name: cond_col: \'Name\', test_col: \'InChI Key\'')
        #     fail_example = _df[~condition_name]
        #     self.logger.debug(self._logging_format_dataframe(fail_example))
        #     clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        _df = clean_df[clean_df['Mixture'] != 'mixture']
        _df = _df.dropna(subset = ['InChI Key'])
        condition_inchikey = _df.apply(lambda x: self._check_inchikey_vs_name(x, cond_col = 'InChI Key', test_col = 'Name', _map = self.map_inchikey_to_synonyms), axis = 1)
        if not condition_inchikey.all():
            passed = False
            self.logger.debug('FAIL in check_inchikey_vs_name: cond_col: \'InChI Key\', test_col: \'Name\'')
            fail_example = _df[~condition_inchikey]
            print(fail_example["InChI Key"].unique())
            print(len(fail_example["InChI Key"].unique()))
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_inchikey_vs_name:  PASS')
        else:
            self.logger.warning('FINISHED: check_inchikey_vs_name:  FAIL')
        return clean_df, passed


    @staticmethod
    def _check_mutation(row, seq_col, mutation_col):
        """
        Check if amino acid that is mutated is in original sequence at the correct position.

        This function is expected to be used with pandas apply.

        Parameters:
        -----------
        row : pandas.Series
            row in the dataframe

        seq_col : str
            name of the column with sequence information

        mutation_col : str
            name of the column containing mutation.

        Returns:
        --------
        passed : bool
            result of the check.
        """
        seq = row[seq_col]
        try:
            mutations = row[mutation_col].strip().replace('\s\s+', '_').split('_')
            passed = True
            for mutation in mutations:
                _from = mutation[0]
                _to = mutation[-1]
                _position = int(mutation[1:-1]) - 1
                if _position >= len(seq):
                    passed = False
                elif not seq[_position] == _from:
                    passed = False
        except Exception as e:
            print(row[mutation_col])
            raise e
        return passed

    def check_mutation(self, full_df):
        """
        Check if amino acid that is supposed to be mutated at a given position can be found on that position.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_mutation')
        passed = True
        clean_df = full_df
        _df = clean_df
        if _df["Mutation"].isna().all():
            passed = True
        else:
            _df = _df.dropna(subset = ['Mutation'])
            condition = _df.apply(lambda x: self._check_mutation(x, seq_col='_Sequence', mutation_col='Mutation'), axis = 1)
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_mutation')
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_mutation:  PASS')
        else:
            self.logger.warning('FINISHED: check_mutation:  FAIL')
        return clean_df, passed


    @staticmethod
    def order_mutations(x, sep):
        if x == x:
            m = x.split(sep)
            m.sort()
            return '_'.join(m)
        else:
            return x

    @staticmethod
    def _check_consistency(x, col):
        return len(x[col].unique()) == 1

    @staticmethod
    def _check_consistency_examples(x, col):
        """
        Get examples that violate check_mutated_sequence_consistency.
        """
        return x[col].unique()

    def check_mutated_sequence_consistency(self, full_df):
        """
        Perform mutation on a sequence and check if we get the same sequence as for some other non-mutated one (or mutated differently).

        Dataframe is first grouped by mutated sequences and we check if there is more than one \"mutated ID\" for them. If yes, check fails.
        Then the dataframe is grouped by \"mutated ID\" and check fails if there is more than 1 sequence for them.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.

        Notes:
        ------
        *This works only for records with UniProt ID.*
        """
        self.logger.info('STARTED: check_mutated_sequence_consistency')
        passed = True
        clean_df = full_df.copy()
        clean_df['mutated_Sequence'] = clean_df.apply(lambda x: perform_mutation(x, mutation_col = 'Mutation', seq_col = '_Sequence'), axis = 1)
        clean_df['ordered_Mutation'] = clean_df['Mutation'].apply(lambda x: self.order_mutations(x, sep = '_'))
        clean_df['mutated_Uniprot ID'] = clean_df.apply(lambda x: x['Uniprot ID'] + '_' + x['ordered_Mutation'] if (isinstance(x['ordered_Mutation'], str)) & (~isinstance(x['Uniprot ID'], float)) else x['Uniprot ID'], axis = 1)#In case of Mutation but no Uniprot ID 
        
        _df = clean_df.copy()
        if not _df['Uniprot ID'].isna().all():
            _df = _df.dropna(subset = ['Uniprot ID'])
        condition_seq = _df.groupby(['mutated_Sequence','species']).apply(lambda x : self._check_consistency(x, col = 'mutated_Uniprot ID'))
        if not condition_seq.all():
            passed = False
            self.logger.debug('FAIL in check_mutated_sequence_consistency: groupby: \'mutated_Sequence\'')
            fail_example = _df[_df['mutated_Sequence'].isin(condition_seq[~condition_seq].index)]
            _fail_example = fail_example.groupby(['mutated_Sequence','species']).apply(lambda x : self._check_consistency_examples(x, col = 'mutated_Uniprot ID'))
            self.logger.debug(self._logging_format_dataframe(fail_example))
            self.logger.debug(_fail_example)
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        _df = clean_df.copy()
        if _df['Uniprot ID'].isna().all():
            passed = True
        else:
            _df = _df.dropna(subset = ['Uniprot ID'])
            condition_id = _df.groupby(['mutated_Uniprot ID']).apply(lambda x : self._check_consistency(x, col = 'mutated_Sequence'))
            if not condition_id.all():
                    passed = False
                    self.logger.debug('FAIL in check_mutated_sequence_consistency: groupby: \'mutated_Uniprot ID\'')
                    fail_example = _df[_df['mutated_Uniprot ID'].isin(condition_id[~condition_id].index)]
                    _fail_example = fail_example.groupby(['mutated_Uniprot ID']).apply(lambda x : self._check_consistency_examples(x, col = 'mutated_Sequence'))
                    self.logger.debug(self._logging_format_dataframe(fail_example))
                    self.logger.debug(_fail_example)
                    clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        
        # TODO: Uncomment if you want to check for the same sequence obtained by different mutations.
        # _df = clean_df.copy()
        # condition_id = _df.groupby('mutated_Sequence').apply(lambda x : self._check_consistency(x, col = 'Mutation'))
        # if not condition_id.all():
        #         passed = False
        #         self.logger.debug('FAIL in check_mutated_sequence_consistency: groupby: \'mutated_Sequence\'')
        #         fail_example = _df[_df['mutated_Sequence'].isin(condition_id[~condition_id].index)]
        #         _fail_example = fail_example.groupby('mutated_Sequence').apply(lambda x : self._check_consistency_examples(x, col = 'Mutation'))
        #         self.logger.debug(self._logging_format_dataframe(fail_example))
        #         self.logger.debug(_fail_example)
        #         clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_mutated_sequence_consistency:  PASS')
        else:
            self.logger.warning('FINISHED: check_mutated_sequence_consistency:  FAIL')
        return clean_df, passed
        

    def check_response_by_article_consistency(self, full_df):
        """
        Check if in the same article they are not claiming oposite Responsivness.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_response_by_article_consistency')
        passed = True
        clean_df = full_df.copy()
        
        _df = clean_df
        _df = _df[_df['Parameter'] == 'ec50']
        if _df.empty:
            passed = True
        else:
            condition = _df.groupby(['mutated_Sequence', 'InChI Key', 'DOI', 'Value_Screen', 'Tag','Cell_line'], dropna=False).apply(lambda x: self._check_consistency(x, col = 'Responsive'))
            if not condition.all():
                passed = False
                self.logger.debug('FAIL in check_response_by_article_consistency ec50')
                _df = _df.reset_index()
                _df = _df.set_index(['mutated_Sequence', 'InChI Key', 'DOI','Value_Screen'])
                fail_example = _df.loc[(condition[~condition].index)]
                fail_example = fail_example.reset_index()
                fail_example = fail_example.set_index('_row_id')
                fail_example.index.name = None
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        
        clean_df = full_df.copy()
        _df = clean_df
        _df = _df[_df['Parameter'] != 'ec50']
        if _df.empty:
            passed = True
        else:
            condition = _df.groupby(['mutated_Sequence', 'InChI Key', 'DOI', 'Value_Screen', 'Tag','Cell_line'], dropna=False).apply(lambda x: self._check_consistency(x, col = 'Responsive'))
            if not condition.all():
                passed = False
                self.logger.debug('FAIL in check_response_by_article_consistency norm and raw')
                _df = _df.reset_index()
                _df = _df.set_index(['mutated_Sequence', 'InChI Key', 'DOI','Value_Screen'])
                fail_example = _df.loc[(condition[~condition].index)]
                fail_example = fail_example.reset_index()
                fail_example = fail_example.set_index('_row_id')
                fail_example.index.name = None
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_response_by_article_consistency:  PASS')
        else:
            self.logger.warning('FINISHED: check_response_by_article_consistency:  FAIL')
        return clean_df, passed

    @staticmethod
    def _check_sep_canonicalSMILES(x, _map):
        passed = True

        if x["InChI Key"] == x["InChI Key"]:
            for inchikey in x["InChI Key"].split(' '):
                if re.search('\.', _map[inchikey]):
                    passed = False
        elif x["canonicalSMILES"] == x["canonicalSMILES"]:
            if re.search('\.', x["canonicalSMILES"]):
                passed = False
        return passed

    def check_sep_canonicalSMILES(self, full_df):
        """
        Check for \'.\' in canonical SMILES

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_canonicalSMILES')
        passed = True
        clean_df = full_df.copy()
        _df = clean_df
        condition = _df.apply(lambda x: self._check_sep_canonicalSMILES(x, _map = self.map_inchikey_to_canonicalSMILES), axis=1)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        
        if passed:
            self.logger.info('FINISHED: check_canonicalSMILES:  PASS')
        else:
            self.logger.warning('FINISHED: check_canonicalSMILES:  FAIL')
        return clean_df, passed

    def check_mixture_based_on_name(self, full_df):
        """
        Check for pattern that can be found in mixture's name. 
        
        For example: 3(2)-mol_name means there is 2-mol_name and/or 3-mol_name.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_mixture_based_on_name')
        passed = True
        clean_df = full_df

        expression = r"[0-9]\([0-9]\)|[0-9]/[0-9]"

        _df = clean_df[clean_df['Mixture'].str.lower() != 'mixture']
        condition_mixname = ~_df['Name'].str.contains(expression)
        if not condition_mixname.all():
            passed = False
            self.logger.debug('FAIL in check_mixture_based_on_name: \'Name\' contains \'{}\' and not a Mixture'.format(expression))
            fail_example = _df[~condition_mixname]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_mixture_based_on_name:  PASS')
        else:
            self.logger.warning('FINISHED: check_mixture_based_on_name:  FAIL')
        return clean_df, passed


    def check_isomers_based_on_name(self, full_df):
        """
        Check for pattern that can be found in isomer's name. 

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_isomers_based_on_name')
        passed = True
        clean_df = full_df

        expression = r"(?i)\([0-9]{0,1}[e,z,s,r](\,[0-9]{0,1}[e,z,s,r]){0,1}\)|cis|trans|^d-|^l-|\(\+\)|\(-\)"

        _df = clean_df[clean_df['Mixture'].str.lower() != 'mixture']
        condition = ~(_df['Name'].str.contains(expression) & _df["InChI Key"].str.contains("UHFFFAOYSA"))
        if not condition.all():
            passed = False
            self.logger.debug('FAIL in check_isomers_based_on_name: \'Name\' contains \'E,Z,R,S,cis,trans,+,-,d or l\' and \'InChI Key\' contains \'UHFFAOYSA\'')
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_isomers_based_on_name:  PASS')
        else:
            self.logger.warning('FINISHED: check_isomers_based_on_name:  FAIL')
        return clean_df, passed


    def check_len_seq(self, full_df):
        """
        Check if sequence is at least 200 amino acids.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_length_sequence')
        passed = True
        clean_df = full_df
        _df = clean_df
        condition = _df["_Sequence"].apply(lambda x: len(x) >= 200 and len(x) <= 380)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug('FAIL in check_length_sequence')
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_length_sequence:  PASS')
        else:
            self.logger.warning('FINISHED: check_length_sequence:  FAIL')
        return clean_df, passed

    def check_ec50_non_zero(self, full_df):
        """
        TODO
        """
        self.logger.info('STARTED: check_ec50_non_zero')
        passed = True
        clean_df = full_df
        _df = clean_df
        _df = _df[_df['Parameter'] == 'ec50']
        if _df.empty:
            passed = True
        else:
            _df = _df[_df['Value'].apply(lambda x: self._check_castable(x, float))]
            condition = _df['Value'].astype(float) != 0.0
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_ec50_non_zero')
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_ec50_non_zero:  PASS')
        else:
            self.logger.warning('FINISHED: check_ec50_non_zero:  FAIL')
        return clean_df, passed


    def check_value_categorical(self, full_df):
        """
        TODO
        """
        self.logger.info('STARTED: check_value_categorical')
        passed = True
        clean_df = full_df
        for case in self.categorical_values:
            _df = clean_df
            if case['col'] == 'Parameter':
                _df = _df.dropna(subset=[case['col']])
            else:
                pass
            _df = _df.dropna(subset=[case['col']])
            condition = _df[case['col']].isin(case['allowed_values'])
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug('FAIL in check_value_categorical:  col: \'{}\'  allowed values: {} '.format(*case.values()))
                self.logger.debug(self._logging_format_dataframe(fail_example))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_value_categorical:  PASS')
        else:
            self.logger.warning('FINISHED: check_value_categorical:  FAIL')
        return clean_df, passed

    def check_mutation_based_on_geneid(self, full_df):
        """
        TODO
        """
        self.logger.info('STARTED: check_mutation_based_on_geneid')
        passed = True
        clean_df = full_df
        _df = clean_df
        _df = _df[(_df["Gene ID"].str.contains(r'_[A-Z]{1}[0-9]+')) & (_df['Sequence'].isna())]
        condition = _df['Mutation'].isna()
        if condition.empty:
            passed = True
        elif condition.all():
            passed = False
            fail_example = _df[condition]
            self.logger.debug('FAIL in check_mutation_based_on_geneid')
            self.logger.debug(self._logging_format_dataframe(fail_example))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_mutation_based_on_geneid:  PASS')
        else:
            self.logger.warning('FINISHED: check_mutation_based_on_geneid:  FAIL')
        return clean_df, passed

    def __call__(self, df):
        """
        Parameters:
        -----------
        df : pandas.DataFrame
            dataframe to process.
        """
        df_join = df.copy()
        df_join.index.name = '_row_id'
        df_join = self.add_implied_columns(df_join)
        df_join, self.check_results['not_nan'] = self.check_not_nan(df_join)
        df_join, self.check_results['check_sep_canonicalSMILES'] = self.check_sep_canonicalSMILES(df_join)
        df_join, self.check_results['conditioned_not_nan'] = self.check_conditioned_not_nan(df_join)
        df_join, self.check_results['castable'] = self.check_castable(df_join)
        df_join, self.check_results['not_castable'] = self.check_not_castable(df_join)
        df_join, self.check_results['format'] = self.check_format(df_join)
        #df_join, self.check_results['mixture_format'] = self.check_mixture_format(df_join)
        df_join, self.check_results['inchikey_on_pubchem'] = self.check_inchikey_on_pubchem(df_join)
        #df_join, self.check_results['check_chirality'] = self.check_chirality(df_join)
        #df_join, self.check_results['inchikey_vs_name'] = self.check_inchikey_vs_name(df_join)
        df_join, self.check_results['mutation'] = self.check_mutation(df_join)
        df_join, self.check_results['mutated_sequence_consistency'] = self.check_mutated_sequence_consistency(df_join)
        #df_join, self.check_results['response_by_article_consistency'] = self.check_response_by_article_consistency(df_join)
        #df_join, self.check_results['mixture_based_on_name'] = self.check_mixture_based_on_name(df_join)
        #df_join, self.check_results['isomers_based_on_name'] = self.check_isomers_based_on_name(df_join)
        df_join, self.check_results['length_sequence'] = self.check_len_seq(df_join)
        df_join, self.check_results['Ec50_non_zero'] = self.check_ec50_non_zero(df_join)
        df_join, self.check_results['value_categorical'] = self.check_value_categorical(df_join)
        df_join, self.check_results['mutation_based_on_geneid'] = self.check_mutation_based_on_geneid(df_join)

        if all(self.check_results.values()):
            self.logger.info('--FINAL STATUS--:\tPASS\n----------------------------\n')
            return True
        else:
            self.logger.critical('--FINAL STATUS--:\tFAIL\n----------------------------\n')
            raise ValueError('--FINAL STATUS--:\tFAIL.  See {} for details'.format(_logging_file_path))



class PostChecker(Checker):
    """
    Perform additional checks after running PostFormatter.

    See documentation of Checker for details about checks.
    """
    def __init__(self, auxillary_dir = 'Data', log_dir = 'logs'):
        self.auxillary_dir = auxillary_dir
        self.log_dir = log_dir

        self._init_logger(__class__.__name__)
        self._init_config()
        # Chage logging columns:
        self.logging_cols.append('Mixture')

        self.check_results = {}

    def __call__(self, df):
        """
        Parameters:
        -----------
        df : pandas.DataFrame
            dataframe to process.
        """
        df_join = df.copy()
        df_join.index.name = '_row_id'
        df_join = self.add_implied_columns(df_join)
        df_join, self.check_results['not_nan'] = self.check_not_nan(df_join)
        df_join, self.check_results['check_sep_canonicalSMILES'] = self.check_sep_canonicalSMILES(df_join)
        df_join, self.check_results['conditioned_not_nan'] = self.check_conditioned_not_nan(df_join)
        df_join, self.check_results['castable'] = self.check_castable(df_join)
        df_join, self.check_results['not_castable'] = self.check_not_castable(df_join)
        df_join, self.check_results['format'] = self.check_format(df_join)
        df_join, self.check_results['mixture_format'] = self.check_mixture_format(df_join)
        df_join, self.check_results['inchikey_on_pubchem'] = self.check_inchikey_on_pubchem(df_join)
        df_join, self.check_results['check_chirality'] = self.check_chirality(df_join)
        #df_join, self.check_results['inchikey_vs_name'] = self.check_inchikey_vs_name(df_join)
        df_join, self.check_results['mutation'] = self.check_mutation(df_join)
        df_join, self.check_results['mutated_sequence_consistency'] = self.check_mutated_sequence_consistency(df_join)
        df_join, self.check_results['response_by_article_consistency'] = self.check_response_by_article_consistency(df_join)
        df_join, self.check_results['mixture_based_on_name'] = self.check_mixture_based_on_name(df_join)
        #df_join, self.check_results['isomers_based_on_name'] = self.check_isomers_based_on_name(df_join)
        df_join, self.check_results['length_sequence'] = self.check_len_seq(df_join)
        df_join, self.check_results['Ec50_non_zero'] = self.check_ec50_non_zero(df_join)
        df_join, self.check_results['value_categorical'] = self.check_value_categorical(df_join)
        df_join, self.check_results['mutation_based_on_geneid'] = self.check_mutation_based_on_geneid(df_join)
        
        if all(self.check_results.values()):
            self.logger.info('--FINAL STATUS--:\tPASS\n----------------------------\n')
            return True
        else:
            self.logger.critical('--FINAL STATUS--:\tFAIL\n----------------------------\n')
            raise ValueError('--FINAL STATUS--:\tFAIL.  See {} for details'.format(_logging_file_path))



class OptionalChecker(Checker):
    """
    Perform optional non mandatory checks, to use after PostFormatter.

    See documentation of Checker for details about checks.
    """
    def __init__(self, auxillary_dir = 'Data', log_dir = 'logs'):
        self.auxillary_dir = auxillary_dir
        self.log_dir = log_dir

        self._init_logger(__class__.__name__)
        self._init_config()
        # Chage logging columns:
        self.logging_cols.append('Mixture')

        self.check_results = {}

    def _load_auxillary(self):
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

        # map_inchikey_to_CID:
        try:
            map_inchikey_to_CID = pandas.read_csv(os.path.join(self.auxillary_dir, 'map_inchikey_to_CID.csv'), sep = ';', index_col = None)
        except FileNotFoundError:
            map_inchikey_to_CID = pandas.DataFrame([], columns = self.map_inchikey_to_CID_cols)
            # raise FileNotFoundError('Support file: \'map_inchikey_to_CID.csv\' not found. Please run pubchem_utils.py to download inchi key to cid mapping.')

        assert len(map_inchikey_to_CID.columns) == len(self.map_inchikey_to_CID_cols)
        for i in range(len(map_inchikey_to_CID.columns)):
            if map_inchikey_to_CID.columns[i] != self.map_inchikey_to_CID_cols[i]:
                raise ValueError('map_inchikey_to_CID has different columns or column positions than in self.map_inchikey_to_CID_cols')
        map_inchikey_to_CID.set_index(self.map_inchikey_to_CID_cols[0], drop = True, inplace = True)
        map_inchikey_to_CID = map_inchikey_to_CID.squeeze()

        # map_inchikey_to_canonicalSMILES:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'r') as jsonfile:
                map_inchikey_to_canonicalSMILES = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_canonicalSMILES = {}

        # map_inchikey_to_synonyms:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_synonyms.json'), 'r') as jsonfile:
                map_inchikey_to_synonyms = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_synonyms = {}

        #map_name_to_inchikeys:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_name_to_inchikeys.json'), 'r') as jsonfile:
                 map_name_to_inchikeys = json.load(jsonfile)
        except FileNotFoundError:
            map_name_to_inchikeys = {}
        
        self.df_uniprot = df_uniprot
        self.map_inchikey_to_CID = map_inchikey_to_CID
        self.map_inchikey_to_canonicalSMILES = map_inchikey_to_canonicalSMILES
        self.map_inchikey_to_synonyms = map_inchikey_to_synonyms
        self.map_name_to_inchikeys = map_name_to_inchikeys
        return df_uniprot, map_inchikey_to_CID, map_inchikey_to_canonicalSMILES, map_inchikey_to_synonyms, map_name_to_inchikeys

    def _update_auxilary_map_name_to_inchikeys(self, full_df):
        map_name_to_inchikeys = self.map_name_to_inchikeys.copy()
        current_idx = pandas.Index(map_name_to_inchikeys.keys(), name = 'Name')
        candidate_idx = full_df[full_df['Mixture'] == 'mono']['Name'].dropna() # TODO: Can we somehow work with mixtures?
        candidate_idx = pandas.Index(candidate_idx.unique())
        new_idx = candidate_idx.difference(current_idx)
        if len(new_idx) > 0:
            self.logger.info('Updating map_name_to_inchikeys...')
            NEW = get_map_name_to_inchikeys(new_idx.tolist())
            map_name_to_inchikeys.update(NEW)
            with open(os.path.join(self.auxillary_dir, 'map_name_to_inchikeys.json'), 'w') as jsonfile:
                json.dump(map_name_to_inchikeys, jsonfile)
        return map_name_to_inchikeys

    
    def _update_auxillary(self, full_df):
        """
        call and other operation necessery for all _update_auxillary_*.
        """
        self.df_uniprot = self._update_auxilary_df_uniprot(full_df)
        self.map_inchikey_to_CID = self._update_auxilary_map_inchikey_to_CID(full_df)
        self.map_inchikey_to_canonicalSMILES = self._update_auxilary_map_inchikey_to_canonicalSMILES(full_df)
        self.map_inchikey_to_synonyms = self._update_auxilary_map_inchikey_to_synonyms(full_df)
        self.map_name_to_inchikey = self._update_auxilary_map_name_to_inchikeys(full_df)
        return 

    def _check_inchikey_vs_name(self, x, cond_col, test_col, _map):
        """
        Check if text in \'test_col\' is in the list of possible synonyms for \'cond_col\'. It is mainly used to check if name of molecule is in synonyms retrieved
        from PubChem by InChI key.

        This is expected to be used with pandas apply.

        Parameters:
        -----------
        x : pandas.Series
            row of a dataframe

        cond_col : str
            name of the column containing identifier. Entry from here are mapped to synonyms using _map.

        test_col : str
            name of the column with text that should be in the list of synonyms.

        _map : dict
            mapping form identifier to synonyms.

        Notes:
        ------
        This works only with non-mixtures only, because mixtures can have name that is unrelated to names of the elements inside.
        *Entries in \'cond_col\' that are not in _map.keys() are ignored*.
        """
        # TODO: This is not checking for cols not in map.keys() ! ! !
        passed = True
        if x[cond_col] in _map.keys():
            x_clean = self._clean_string_name(x[test_col])
            
            for synonym in _map[x[cond_col]]:
                try: 
                    self._clean_string_name(synonym)
                except Exception as e:
                    print(_map[x[cond_col]], x[cond_col])
                    raise e
            if not x_clean in [self._clean_string_name(synonym) for synonym in _map[x[cond_col]]]:
                passed = False
        return passed

    def check_inchikey_vs_name(self, full_df):
        """
        Check if names can be found inside synonyms retrieved by InChI Key.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.

        Notes:
        ------
        This works only with non-mixtures only, because mixtures can have name that is unrelated to names of the elements inside.
        """
        # TODO: For now this is checked only for non-mixtures.
        self.logger.info('STARTED: check_inchikey_vs_name')
        passed = True
        clean_df = full_df.copy()

        _df = clean_df[clean_df['Mixture'] != 'mixture']
        _df = _df.dropna(subset = ['Name'])
        condition_name = _df.apply(lambda x: self._check_inchikey_vs_name(x, cond_col = 'Name', test_col = 'InChI Key', _map = self.map_name_to_inchikeys), axis = 1)
        if not condition_name.all():
            passed = False
            self.logger.debug('FAIL in check_inchikey_vs_name: cond_col: \'Name\', test_col: \'InChI Key\'')
            fail_example = _df[~condition_name]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset=['Name'])))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        _df = clean_df[clean_df['Mixture'] != 'mixture']
        _df = _df.dropna(subset = ['InChI Key'])
        condition_inchikey = _df.apply(lambda x: self._check_inchikey_vs_name(x, cond_col = 'InChI Key', test_col = 'Name', _map = self.map_inchikey_to_synonyms), axis = 1)
        if not condition_inchikey.all():
            passed = False
            self.logger.debug('FAIL in check_inchikey_vs_name: cond_col: \'InChI Key\', test_col: \'Name\'')
            fail_example = _df[~condition_inchikey]
            print(fail_example["InChI Key"].unique())
            print(len(fail_example["InChI Key"].unique()))
            self.logger.debug(self._logging_format_dataframe(fail_example))
            self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset=['Name'])))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_inchikey_vs_name:  PASS')
        else:
            self.logger.warning('FINISHED: check_inchikey_vs_name:  FAIL')
        return full_df, passed

    @staticmethod
    def _check_chirality(x, _map):
        """
        Check if molecules that are supposed to be achiral have more than 1 stereoisomer. 

        InChI keys containing \'-UHFFFAOYSA-\' should be achiral as well as all having cacnonical SMILES.

        Parameters:
        -----------
        x : pandas.Series
            row of the dataframe.

        _map : dict
            mapping from InChI key to canonical SMILES.
        """
        passed = True
        if x['InChI Key'] == x['InChI Key']:
            for inchikey in x['InChI Key'].split(' '):
                if '-UHFFFAOYSA-' in inchikey:
                    isomers = enumerate_isomers(_map[inchikey])
                    # mol = Chem.MolFromSmiles(_map[inchikey])
                    # chiralCenters = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True, includeCIP=False, useLegacyImplementation=False)
                    # isomers = tuple(EnumerateStereoisomers(mol))
                    if len(isomers) > 1:
                        passed = False
        elif x["canonicalSMILES"] == x["canonicalSMILES"]:
            isomers = enumerate_isomers(x["canonicalSMILES"])
            if len(isomers) > 1:
                passed = False
        return passed

    def check_chirality(self, full_df):
        """
        For InChI keys containing \'-UHFFFAOYSA-\' and for canonical SMILES check if number of stereoisomers is 1. These records should not have isomers
        and check fails if there are any (i.e. more than 1). For multiple elements the check is performed for each element and separator is assumed to
        be *space*. NaNs are ignored. Isomers are identified using rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers (see enumerate_isomers is utils.py).

        InChI keys are mapped to canonical SMILES using \'map_inchikey_to_canonicalSMILES\' and stereoisomers are identified using these SMILES.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.

        Notes:
        ------
        Doesn't work on mixtures.
        """
        self.logger.info('STARTED: check_chirality')
        passed = True
        clean_df = full_df.copy()

        _df = clean_df[clean_df['Mixture'] == "mixture"]
        if _df.empty:
            passed = True
        else:
            condition = _df.apply(lambda x: self._check_chirality(x, _map = self.map_inchikey_to_canonicalSMILES), axis=1)
            if not condition.all():
                passed = False
                fail_example = _df[~condition]
                self.logger.debug(self._logging_format_dataframe(fail_example))
                self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset=['Name'])))
                clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]

        if passed:
            self.logger.info('FINISHED: check_chirality:  PASS')
        else:
            self.logger.warning('FINISHED: check_chirality:  FAIL')
        return full_df, passed

    def check_len_seq(self, full_df):
        """
        Check if sequence is at least 300 and less than 330 amino acids.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_length_sequence')
        passed = True
        clean_df = full_df
        _df = clean_df
        condition = _df["_Sequence"].apply(lambda x: len(x) >= 300 and len(x) <= 330)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug('FAIL in check_length_sequence')
            self.logger.debug(self._logging_format_dataframe(fail_example))
            self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset=['Name'])))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_length_sequence:  PASS')
        else:
            self.logger.warning('FINISHED: check_length_sequence:  FAIL')
        return full_df, passed

    def check_isomers_based_on_name(self, full_df):
        """
        Check for pattern that can be found in isomer's name. 

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_isomers_based_on_name')
        passed = True
        clean_df = full_df

        expression = r"(?i)\([0-9]{0,1}[e,z,s,r](\,[0-9]{0,1}[e,z,s,r]){0,1}\)|cis|trans|^d-|^l-|\(\+\)|\(-\)"

        _df = clean_df[clean_df['Mixture'].str.lower() != 'mixture']
        condition = ~(_df['Name'].str.contains(expression) & _df["InChI Key"].str.contains("UHFFFAOYSA"))
        if not condition.all():
            passed = False
            self.logger.debug('FAIL in check_isomers_based_on_name: \'Name\' contains \'E,Z,R,S,cis,trans,+,-,d or l\' and \'InChI Key\' contains \'UHFFAOYSA\'')
            fail_example = _df[~condition]
            self.logger.debug(self._logging_format_dataframe(fail_example))
            self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset=['Name'])))
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_isomers_based_on_name:  PASS')
        else:
            self.logger.warning('FINISHED: check_isomers_based_on_name:  FAIL')
        return full_df, passed

    def check_blast_result(self, full_df):
        """
        Check if sequence identity to blast sequence is at least 96%.

        Parameters:
        -----------
        full_df : pandas.DataFrame
            dataframe being processed

        Returns:
        --------
        clean_df : pandas.DataFrame
            dataframe where rows that do not pass the check are discarded.

        passed : bool
            result of the check.
        """
        self.logger.info('STARTED: check_blast_result')
        passed = True
        clean_df = full_df
        _df = clean_df
        condition = _df["blast_identity"].apply(lambda x: x >= 96)
        if not condition.all():
            passed = False
            fail_example = _df[~condition]
            self.logger.debug('FAIL in check_blast_result')
            self.logger.debug(self._logging_format_dataframe(fail_example.drop_duplicates(subset='mutated_Sequence')))
            self.logger.debug(fail_example.drop_duplicates(subset='mutated_Sequence')[['Gene ID','Uniprot ID','blast_uniprot_id','species','mutated_Sequence','blast_seq']])
            clean_df = clean_df.loc[clean_df.index.difference(fail_example.index)]
        if passed:
            self.logger.info('FINISHED: check_blast_result:  PASS')
        else:
            self.logger.warning('FINISHED: check_blast_result:  FAIL')
            return full_df, passed

    def __call__(self, df):
        """
        Parameters:
        -----------
        df : pandas.DataFrame
            dataframe to process.
        """
        df_join = df.copy()
        df_join.index.name = '_row_id'
        df_join = self.add_implied_columns(df_join)
        df_join, self.check_results['check_blast_result'] = self.check_blast_result(df_join)
        df_join, self.check_results['check_chirality'] = self.check_chirality(df_join)
        df_join, self.check_results['inchikey_vs_name'] = self.check_inchikey_vs_name(df_join)
        df_join, self.check_results['length_sequence'] = self.check_len_seq(df_join)
        df_join, self.check_results['isomers_based_on_name'] = self.check_isomers_based_on_name(df_join)
        
        if all(self.check_results.values()):
            self.logger.info('--FINAL STATUS--:\tPASS\n----------------------------\n')
            return True
        else:
            self.logger.critical('--FINAL STATUS--:\tFAIL\n----------------------------\n')
            raise ValueError('--FINAL STATUS--:\tFAIL.  See {} for details'.format(_logging_file_path))



if __name__ == '__main__':
    from formatting import PreFormatter
    formatter = PreFormatter()
    checker = Checker()

    df = pandas.read_csv("./RawData/Chemodb_20220311.csv", sep=';', index_col = None)

    # print(df.columns)
    df, exclude_df = formatter(df)
    checker(df)

    # print(df.loc[8039]['InChI Key'])
    # print(df.loc[8039]['InChI Key'].replace(r'  ', ''))

    # for s in df.loc[8039]['InChI Key']:
    #     print('{} : {}'.format(s, ord(s)))