import os
import sys
import pandas
import json
import itertools
import logging

from utils import enumerate_isomers
from pubchem_utils import get_map_inchikey_to_canonicalSMILES, get_map_isomericSMILES_to_inchikey



class IsoRetriever:
    def __init__(self):
        self.auxillary_dir = 'Data'
        self.log_dir = 'logs'

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

        self._load_auxillary()


    def _load_auxillary(self):
        # map_inchikey_to_canonicalSMILES:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_inchikey_to_canonicalSMILES.json'), 'r') as jsonfile:
                map_inchikey_to_canonicalSMILES = json.load(jsonfile)
        except FileNotFoundError:
            map_inchikey_to_canonicalSMILES = {}

        # map_isomericSMILES_to_inchikey:
        try:
            with open(os.path.join(self.auxillary_dir, 'map_isomericSMILES_to_inchikey.json'), 'r') as jsonfile:
                map_isomericSMILES_to_inchikey = json.load(jsonfile)
        except FileNotFoundError:
            map_isomericSMILES_to_inchikey = {}

        #
        self.map_inchikey_to_canonicalSMILES = map_inchikey_to_canonicalSMILES
        self.map_isomericSMILES_to_inchikey = map_isomericSMILES_to_inchikey
        return map_inchikey_to_canonicalSMILES, map_isomericSMILES_to_inchikey


    def _update_auxilary_map_inchikey_to_canonincalSMILES(self, candidate_idx):
        print(candidate_idx)
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
        self.map_inchikey_to_canonicalSMILES = map_inchikey_to_canonicalSMILES
        return map_inchikey_to_canonicalSMILES


    def _update_auxilary_map_isomericSMILES_to_inchikey(self, candidate_idx):
        assert candidate_idx.name == 'isomericSMILES'
        map_isomericSMILES_to_inchikey = self.map_isomericSMILES_to_inchikey.copy()
        current_idx = pandas.Index(map_isomericSMILES_to_inchikey.keys(), name = 'isomericSMILES')
        new_idx = candidate_idx.difference(current_idx)
        if len(new_idx) > 0:
            self.logger.info('Updating map_isomericSMILES_to_inchikey...')
            NEW = get_map_isomericSMILES_to_inchikey(new_idx.tolist())
            map_isomericSMILES_to_inchikey.update(NEW)
            with open(os.path.join(self.auxillary_dir, 'map_isomericSMILES_to_inchikey.json'), 'w') as jsonfile:
                json.dump(map_isomericSMILES_to_inchikey, jsonfile)
        self.map_isomericSMILES_to_inchikey = map_isomericSMILES_to_inchikey
        return map_isomericSMILES_to_inchikey


    def retrieve_isomericSMILES(self, df):
        _df = df.copy()
        _df = _df.dropna(subset = ['InChI Key'])
        candidate_idx = _df['InChI Key'].dropna()
        candidate_idx = candidate_idx.str.split(' ').explode() # TODO: pandas FutureWarning for this row.
        candidate_idx = pandas.Index(candidate_idx.unique(), name = 'InChI Key')
        self._update_auxilary_map_inchikey_to_canonincalSMILES(candidate_idx)
        map_inchikey_to_isomers = {}
        for x in _df['InChI Key'].unique():
            _res = {}
            for inchikey in x.split(' '):
                if '-UHFFFAOYSA-' in inchikey: # This is done only for non-isomeric
                    _res[inchikey] = enumerate_isomers(self.map_inchikey_to_canonicalSMILES[inchikey])
            map_inchikey_to_isomers.update(_res)

        isomericSMILES = list(itertools.chain(*map_inchikey_to_isomers.values()))
        candidate_idx = pandas.Index(isomericSMILES, name = 'isomericSMILES')
        candidate_idx = candidate_idx.unique()
        self._update_auxilary_map_isomericSMILES_to_inchikey(candidate_idx)
        self.map_inchikey_to_isomers = map_inchikey_to_isomers
        return self.map_isomericSMILES_to_inchikey


    def retrieve_isoInChIKey(self, df):
        self.retrieve_isomericSMILES(df)
        map_inchikey_to_isomers = self.map_inchikey_to_isomers.copy()
        map_inchikey_to_isoinchikey = {key : list(map(self.map_isomericSMILES_to_inchikey.get, val)) for key, val in map_inchikey_to_isomers.items()}
        self.map_inchikey_to_isoinchikey = map_inchikey_to_isoinchikey
        return self.map_inchikey_to_isoinchikey

    # TODO: is canonicalSMILES treated as well? 


if __name__ == '__main__':
    from formatting import PreFormatter
    df = pandas.read_csv('/home/matej/Desktop/tmp/Chemodb_20210110.csv', sep=';', index_col = 0)

    formatter = PreFormatter()
    test, exclude_test = formatter(df)

    post_formatter = IsoRetriever()
    post_formatter.retrieve_isoInChIKey(test)