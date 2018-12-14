# Physical Sciences in Oncology cell line characterization dataset

## Introduction

The cell-line characterization project (see https://physics.cancer.gov/ and https://physics.cancer.gov/data/ for more information on PS-ON and related data resources) integrates cell physical-phenotypes with omics data, capturing all available relationships between experimental measurements of cell physical attributes (e.g. cell area, volume, circularity, nuclear volume, stiffness, contractility, etc.) and various omics data (e.g. exome, mRNA, and miRNA sequencing and proteomics assays), across cell-culture substrates and cell lines. For more information on available datasets, raw data files, and experimental measurements, please refer to the data center coordination portal at https://nciphub.org/groups/nci_physci/psondcc. 

Additional resources on related studies, experiments as well as curated raw and processed data are accessible at https://www.synapse.org/#!Synapse:syn7248578/wiki/405995.

## Data model

A completely specified data model describing the available cell-line characterization data is available [here](https://github.com/Sage-Bionetworks/PSON-cell-line-panel-characterization/blob/master/db-scheme.pdf). See figures A, B and C below for a brief data-model overview below. 

![alt text](https://raw.githubusercontent.com/milen-sage/PSON-cell-line-panel-characterization/master/descriptive_data_model_fig.png)

Physical phenotype measurements can be looked up across all (or a subset) of cell lines, all (or a subset) of culture substrates, and can be cross-referenced with any available omics measurements (e.g. gene expression) in the corresponding cell lines and culture substrates. The data model is readily instantiable as a relational database (RDB), but supports other implementations as well (e.g. stand-alone data frames - https://www.synapse.org/#!Synapse:syn17096709); please refer to the following section for more details on available implementations, interfaces, and data model use. 

The cell-line panel characterization data is integrated within a relation-oriented data model. A. Physical cell measurements comprise five categories of relational entities; omics data comprise three categories. The data model captures all available relationships (solid lines) between physical phenotype measurements (left), omics measurements (right), the culture substrates and cell lines for which measurements were performed (center). B. Top: each relational entity in the data model (examples highlighted in orange lines), include single-cell measurements (or summary statistics of single-cell measurements) across different cell attributes related to a given phenotype category (e.g. motility, proteomics, etc.). All relational entities (and their respective attributes) containing measurements data are linked to the corresponding cell line and culture substrate attributes for which the measurements were made. The cell line and substrate relational entities (green) provide further characteristics of the experiment environment (e.g. culture gel type, stiffness, etc.) and the reagent cell line (e.g. source tissue, diagnosis, etc.). E.g., motility attributes can be cross-linked to proteomics data for all (or a subset of) gel types (e.g. Hyaloronic acid matrix) and cell lines (e.g. A375, skin tissue cancer) in which measurements were performed, which facilitates correlation analysis. Bottom: each relational entity corresponds to a relational table following a provided database schema implementing the data model. A relational database implementation (along with dataframe-based implementation) of the data model is available in an open-access repository [see links]. Various programming interfaces can be used to query the data, joining on all (or subset) of available facets.

## Data interfaces

All data can be accessed at https://www.synapse.org/#!Synapse:syn7248578/wiki/405995. Data can be filtered across various facets to obtain relevant experimental measurements and analytic results.

In addition, integrated datasets following the above data model are provided in the form of 

1) data frames - https://www.synapse.org/#!Synapse:syn17096709
2) query-able data tables - https://www.synapse.org/#!Synapse:syn17096706/tables/
3) SQLite database - https://www.synapse.org/#!Synapse:syn17096714
4) MySQL database can be readily instantiated (with included field annotations and comments):
    - download data from Synapse (free account creation is required and takes 3 minutes): https://www.synapse.org/#!Synapse:syn17096709
    - clone this repository
    - assuming a working MySQL installation, login to mysql and create a new database

        mysql -u username -p 

        mysql>CREATE DATABASE "dbname";

        mysql>USE "dbname";

    - run the SQL script db-scheme-sql.sql

        mysql>source path/to/db-scheme-sql.sql

    - run the SQL script db-populate.sql (change data source path in the script to point to the downloaded dataframes)
    
        mysql>source path/to/db-populate.sql

## Data usage

All data interfaces follow the provided database [schema](https://github.com/Sage-Bionetworks/PSON-cell-line-panel-characterization/blob/master/db-scheme.pdf).

These instructions are continuously updated. Sample queries in Python, R and SQL are provided below (note, only relevant queries are provided; these are not fully functional scripts as data needs to be ingested):

R: library(dplyr); left_join(motility, proteomics_iBAQ, by = c(‘cl_id’, ‘sub_id’));

Python: import pandas; motility.merge(proteomics_iBAQ, how = ‘left’, on =[‘cl_id’, ‘sub_id’]) 

SQL: SELECT * FROM motility LEFT JOIN proteomics_iBAQ ON motility.cl_id = proteomics_iBAQ.cl_id AND motility.sub_id = proteomics_iBAQ.sub_id;
