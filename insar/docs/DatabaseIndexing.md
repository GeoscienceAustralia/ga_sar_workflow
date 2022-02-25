## InSAR geospatial database

For large scale datasets, it's not feasible to individually pass in data acquisitions and have their data parsed every time the stack needs to start/resume/append processing data.  To this end, `gamma_insar` supports indexing data acquisitions into a geospatial temporal database which keeps a record of all the data acquisitions available, their geospatial metadata, and indexes them in a way that's optimal for the workflow to query and use.

This guide explains how to generate and maintain such a database.

#### Metadata YAML extraction

The first step in creating the database is to extract the appropriate metadata from acquisitions into a standardised format, this must be done for all acquisitions you wish to add to the database.

This example extracts the SLC acquisition details into YAML files for a single month.  It takes 1-2 hours for ~5 years of Sentinel-1 acquisitions.

```BASH
slc-archive slc-ingestion \
    --save-yaml \
    --yaml-dir <dir-to-save-generated-yamls> \
    --year <year-to-process> \
    --month <month-to-process> \
    --slc-dir /g/data/fj7/Copernicus/Sentinel-1/C-SAR/SLC \  # change if required
    --log-pathname <filename-to-save-log-to>
```

#### Database Creation

Once the metadata has been standardised into the yaml format, it can be ingested into the database.

This example creates (or updates an existing) the database (sqlite .db) file from SLC metadata stored in the yaml files:

```BASH
slc-archive slc-ingest-yaml \
        --database-name <output-sqlite-database-filename> \
        --yaml-dir <base-dir-containing-yaml-files> \
        --log-pathname <output-json-log-filename>
```

Note: the run time can be quite long depending on how many yaml files there are.
