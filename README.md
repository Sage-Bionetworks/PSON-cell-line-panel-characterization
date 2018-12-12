# PSON-cell-line-panel-characterization
PSON cell line data model schema, database scripts, and related utilities

To instantiate a MySQL database:

1) download data from Synapse (free account creation is required and takes 5 minutes):
https://www.synapse.org/#!Synapse:syn17093184

2) clone this repository

3) assuming a working MySQL installation, login to mysql on the command line and create a new database

mysql -u root -p
mysql>CREATE DATABASE "dbname";
mysql>USE "dbname";

4) Run the SQL script db-scheme-sql.sql
mysql>source path/to/db-scheme-sql.sql

5) Run the SQL script db-populate.sql
mysql>source path/to/db-populate.sql

To instantiate a SQLite database:

1) download data from Synapse (free account creation is required and takes 5 minutes):
https://www.synapse.org/#!Synapse:syn17093184

2) clone this repository

3) run create_sql_db.py on command line
python create_sql_db.py