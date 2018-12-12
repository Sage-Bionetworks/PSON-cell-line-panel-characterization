import os
import pandas as pd
import synapseclient
from synapseclient import File
import glob
import sqlite3 as sq
from sqlalchemy import create_engine

syn = synapseclient.login()

# project temp workspace
db_project = 'syn17093184'

# database csv tables path
db_tables_path = "/Users/milen-sage/workspace/PSON-cell-line-data/schema_tables"
sqlite_store_path = os.path.join(db_tables_path, "pson_cell_line_panel_characterization.sqlite")

table_files = glob.glob(os.path.join(db_tables_path, "*.csv"))

table_names = [os.path.basename(path).split(".")[0] for path in table_files]
tab_separated = ["miRNA_mapping_annotation.csv", "miRNA_mapping_summary.csv", "mRNA_gene_transcript_annotation.csv", "proteomics_iBAQ.csv", "proteomics_phospho.csv", "proteomics_tmt.csv"]


con = sq.connect(sqlite_store_path)
#cur = con.cursor()

for table_file in table_files:
    table_name = os.path.basename(table_file).split(".")[0] 
    print(table_name)
    if os.path.basename(table_file) in tab_separated:
        df = pd.read_csv(table_file, sep = '\t')
    else:
        df = pd.read_csv(table_file, sep = ',')
    df.rename(columns = {"condition_id":"sub_id", "condition":"substrate"}, inplace = True)
    print("Processed file")

    print("Creating table...")

    df.to_sql(table_name, con, if_exists='replace', index = False) 
    print("Done.")

con.commit()
con.close()


sqlite_db_file = File(sqlite_store_path, description='PSON cell line panel characterization SQLite DB', parent = db_project)
sqlite_db_file = syn.store(sqlite_db_file)
