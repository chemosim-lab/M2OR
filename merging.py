import os
import pandas
import psycopg
import json

from psycopg.rows import dict_row

from configs import config_db

# (OK) TODO: Create sequence_id - This is unnecessary
# (OK) TODO: get CAS ... do we need it?? - yes, but only as user input and we will not search for it.
# TODO: get CID

class MergerSQL:
    """
    class to work with PostgreSQL database. It is used for merging new data and exporting from database.
    
    The main purpose of using SQL is to have nontrivial UNIQUE condition on columns.

    Attributes:
    -----------
    config_file : str
        path to Postgre init file.

    conn : psycopg.Connection
        connection to database.

    cur : psycopg.Cursor
        database cursor.

    Notes:
    ------
    see \'psql\' command in linux.
    """
    def __init__(self, config_file):
        """
        Parameters:
        -----------
        config_file : str
            path to Postgre init file.
        """
        self.config_file = config_file


    def _connect_to_db(self):
        """
        Create connection and cursor for database.

        References:
        -----------
        (some tutorials) \n
        https://www.postgresqltutorial.com/postgresql-python/connect/ \n
        https://www.psycopg.org/psycopg3/docs/basic/usage.html \n
        https://pythontic.com/database/postgresql/create%20database \n
        """
        params = config_db(filename = self.config_file)
        self.conn = psycopg.connect(**params, row_factory=dict_row)
        self.cur = self.conn.cursor()
        return None


    def open(self):
        """
        Create connection and cursor for database.
        """
        self._connect_to_db()


    def create_table(self):
        """
        Create SQL table with name \'main\' and add UNIQUE condition to it. 
        See code for column names and types.

        Empty columns (i.e. None) is considered as equal in UNIQUE condition and the combination of
        ["Mutation", "Uniprot ID", "Sequence", "InChI Key", "canonicalSMILES",  "Parameter",  
        "Value",  "Unit",  "Value_Screen",  "Unit_Screen",  "Responsive",  "Type", "Tag",  "DOI"] 
        is considered unique.
        """
        # TODO: This is here for debug only. Delete afterwards.
        # self.cur.execute("""
        #     DROP TABLE IF EXISTS main
        #     """)
        # self.cur.execute("""
        #     DROP TABLE main_test
        #     """)

        self.cur.execute("""
            CREATE TABLE main (
                id                  serial PRIMARY KEY,
                "species"             text,
                "Mutation"            text,
                "Gene ID"             text,
                "Uniprot ID"          varchar(30),
                "Sequence"            text,
                "Name"                text,
                "CID"                 text,
                "CAS"                 text,
                "InChI Key"           text,
                "canonicalSMILES"     text,
                "Parameter"           text,
                "Value"               text,
                "Unit"                varchar(10),
                "Value_Screen"        text,
                "Unit_Screen"         varchar(10),
                "Responsive"          integer,
                "nbr_measurements"    integer,
                "Type"                text,
                "Cell_line"           text,
                "Delivery"            text,
                "Assay"               text,
                "Gprotein"            text,
                "Co_transfection"     text,
                "Assay System"        text,
                "Tag"                 text,
                "Reference"           text,
                "DOI"                 text,
                "Reference Position"  text,
                "Mixture"             varchar(30),
                "Norm_Foreign_Key"    text, 
                CHECK("Mutation" <> 'X9999X'),
                CHECK("Uniprot ID" <> '0000000000'),
                CHECK("Sequence" <> 'XXXXX'),
                CHECK("InChI Key" <> 'XXXXXXXXXXXXXX-XXXXXXXXXX-X'),
                CHECK("canonicalSMILES" <> '0000000000'),
                CHECK("Value" <> '0000000000'),
                CHECK("Unit" <> '0000000000'),
                CHECK("Value_Screen" <> '0000000000'),
                CHECK("Unit_Screen" <> '0000000000'),
                CHECK("Type" <> '0000000000'),
                CHECK("Tag" <> '0000000000')
                )
            """)
            # "seq_id"              text NOT NULL,

        self.cur.execute("""
                    CREATE UNIQUE INDEX 
                        main_unique_entry_index ON main
                    (COALESCE("Mutation", 'X9999X'),
                    COALESCE("Uniprot ID", '0000000000'),
                    COALESCE("Sequence", 'XXXXX'),
                    COALESCE("InChI Key", 'XXXXXXXXXXXXXX-XXXXXXXXXX-X'),
                    COALESCE("canonicalSMILES", '0000000000'), 
                    "Parameter", 
                    COALESCE("Value", '0000000000'), 
                    COALESCE("Unit", '0000000000'), 
                    COALESCE("Value_Screen", '0000000000'), 
                    COALESCE("Unit_Screen", '0000000000'), 
                    "Responsive", 
                    COALESCE("Type", '0000000000'),
                    COALESCE("Tag", '0000000000'), 
                    "DOI");
        """)

        self.conn.commit()
        return None


    def update_db(self, new_df):
        """
        Insert new records to the database. If the record already exists (i.e. violation of UNIQUE condition)
        the old record is kept and new one is ignored (i.e. \'ON CONFLICT DO NOTHING\' in SQL).
        """
        df = new_df.copy()
        df = df.replace({float('nan') : None})
        rows = df.to_dict(orient = 'records')

        self.cur.executemany("""
            INSERT INTO main ("species", "Mutation", "Gene ID", "Uniprot ID", "Sequence", 
                            "Name", "CID", "CAS", "InChI Key", "canonicalSMILES", "Parameter", "Value", "Unit", 
                            "Value_Screen", "Unit_Screen", "Responsive", "nbr_measurements", "Type", 
                            "Cell_line", "Delivery", "Assay", "Gprotein", "Co_transfection", "Assay System", "Tag",
                            "Reference", "DOI", "Reference Position", "Mixture") 
                        VALUES (%(species)s, %(Mutation)s, %(Gene ID)s, %(Uniprot ID)s, %(Sequence)s, 
                                %(Name)s, %(CID)s, %(CAS)s, %(InChI Key)s, %(canonicalSMILES)s, %(Parameter)s, %(Value)s, %(Unit)s, 
                                %(Value_Screen)s, %(Unit_Screen)s, %(Responsive)s, %(nbr_measurements)s, %(Type)s, 
                                %(Cell_line)s, %(Delivery)s, %(Assay)s, %(Gprotein)s, %(Co_transfection)s, %(Assay System)s, %(Tag)s,
                                %(Reference)s, %(DOI)s, %(Reference Position)s, %(Mixture)s) ON CONFLICT DO NOTHING
        """, rows)

        self.conn.commit()


    def DEBUG_update_db(self, new_df):
        """
        TODO: For DEBUG only. This prints out some rows with violation of UNIQUE (at least I hope it does). Delete afterwards...
        """
        df = new_df.copy()
        df = df.replace({float('nan') : None})
        rows = df.to_dict(orient = 'records')

        i = 0
        for row in rows:
            try:
                self.cur.execute("""
                INSERT INTO main ("species", "Mutation", "Gene ID", "Uniprot ID", "Sequence", 
                            "Name", "CID", "CAS", "InChI Key", "canonicalSMILES", "Parameter", "Value", "Unit", 
                            "Value_Screen", "Unit_Screen", "Responsive", "nbr_measurements", "Type", 
                            "Cell_line", "Co_transfection", "Assay System", "Tag",
                            "Reference", "DOI", "Reference Position", "Mixture") 
                        VALUES (%(species)s, %(Mutation)s, %(Gene ID)s, %(Uniprot ID)s, %(Sequence)s, 
                                %(Name)s, %(CID)s, %(CAS)s, %(InChI Key)s, %(canonicalSMILES)s, %(Parameter)s, %(Value)s, %(Unit)s, 
                                %(Value_Screen)s, %(Unit_Screen)s, %(Responsive)s, %(nbr_measurements)s, %(Type)s, 
                                %(Cell_line)s, %(Co_transfection)s, %(Assay System)s, %(Tag)s,
                                %(Reference)s, %(DOI)s, %(Reference Position)s, %(Mixture)s) 
                """, row) # ON CONFLICT DO NOTHING

            except Exception as e:
                print(row['index'])
                print(row['Sequence'])
                print('--------')
                print(i)
                if i > 2:
                    raise e
                i += 1

        
    def export(self):
        """
        Export entire \'main\' table to csv separated by \';\'.
        """
        self.cur.execute("SELECT * FROM main")
        df = self.cur.fetchall()
        df = pandas.DataFrame.from_records(df)
        # df = df.replace({None : float('nan')})
        return df

    def close(self):
        """
        Close database cursor and connection.
        """
        self.cur.close()
        self.conn.close()




if __name__ == '__main__':
    
    from formatting import PreFormatter,PostFormatter
    from checking import Checker,PostChecker

    df = pandas.read_csv('chemodb_20230320_tocheck.csv', sep=';', index_col = None)
    
    formatter = PreFormatter()
    df, exclude_df = formatter(df)
    # exclude_df.to_csv('RawData_test/exclude.csv', sep=';')
    checker = Checker()
    checker(df)
    postformatter = PostFormatter()
    df = postformatter(df)
    postchecker = PostChecker()
    postchecker(df)
    config_file = 'Database/chemosimdb.ini'
    merger = MergerSQL(config_file)
    merger.open()
    merger.create_table()
    merger.DEBUG_update_db(df)
    df = merger.export()
    merger.close()
    # [9519 rows x 21 columns]

    print(df)

