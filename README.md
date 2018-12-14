# Physical Sciences in Oncology cell line characterization dataset

## Introduction

The cell-line characterization project (see https://physics.cancer.gov/ and https://physics.cancer.gov/data/ for more information on PS-ON and related data resources) integrates cell physical-phenotypes with omics data, capturing all available relationships between experimental measurements of cell physical attributes (e.g. cell area, volume, circularity, nuclear volume, stiffness, contractility, etc.) and various omics data (e.g. exome, mRNA, and miRNA sequencing and proteomics assays), across cell-culture substrates and cell lines. For more information on available datasets, raw data files, and experimental measurements, please refer to the data center coordination portal at https://nciphub.org/groups/nci_physci/psondcc. 

Additional resources on related studies, experiments as well as curated raw and processed data are accessible at https://www.synapse.org/#!Synapse:syn7248578/wiki/405995.

## Data model

A completely specified data model describing the available cell-line characterization data is available [here](https://github.com/Sage-Bionetworks/PSON-cell-line-panel-characterization/blob/master/db-scheme.pdf). See figures A, B and C below for a brief data-model overview below. 

![alt text](https://raw.githubusercontent.com/milen-sage/PSON-cell-line-panel-characterization/master/descriptive_data_model_fig.png)

Physical phenotype measurements can be looked up across all (or a subset) of cell lines, all (or a subset) of culture substrates, and can be cross-referenced with any available omics measurements (e.g. gene expression) in the corresponding cell lines and culture substrates. The data model is readily instantiable as a relational database (RDB), but supports other implementations as well (e.g. stand-alone data frames - https://www.synapse.org/#!Synapse:syn17096709); please refer to the following section for more details on available implementations, interfaces, and data model use. 

The data model organizes all experimental data in five categories of physical phenotype properties and three categories of omics data (figure A). Phenotype categories deriving from single-cell physical measurements include cell motility, stiffness, morphology, contractility, and proliferation. Each of these categories encompasses a set of cell attributes (e.g. cell motility is characterized by cell speed, position, etc.) and corresponds to a single relational entity (e.g. a RDB table). The relational entity stores all single-cell measurements for each cell attribute within the corresponding category. As applicable, a separate relational entity provides summary statistics over all single-cell measurements for a physical attribute (e.g. a separate RDB table for mean and standard deviation of motility measurements). Each entry in a relational entity (e.g. cell speed measurement) is linked to the culture substrate and cell line used for the corresponding experiment. Analogously, omics data categories deriving from various sequencing and proteomics assays included genomics, proteomics, and RNA measurements. Within an omics category, measurements stored in each relation entity (e.g. iBAQ proteomics) are also linked to the culture substrate and cell line used for the corresponding assay. Each relational entity is annotated using standardardized controlled vocabulary terms and metadata, as available. For instance, variant calls annotation follows mutation analysis file (MAF) Genomic Data Commons (GDC) terms - allowing interfacing with other existing data resources. 

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
