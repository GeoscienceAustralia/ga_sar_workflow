## GAMMA-INSAR

A tool to process Sentinel-1 SLC to Analysis Ready Data (ARD) using GAMMA SOFTWARE. The ARD products are georeferenced backscatter and interferograms.

Using `gamma_insar` currentlty requires a GitHub account (& GA access permissions as it's not yet open source), as we do not currently deploy
releases any other way than via github.

#### Operating System tested

* Linux
  * CentOS Linux release 8.3.2011 (NCI's Gadi environment)
  * Ubuntu 18.04 flavours (unit tests only)
  * [osgeo/gdal docker images](https://hub.docker.com/r/osgeo/gdal) (unit tests only)
* macOS (unit tests only)
  * macOS v11.x (GDAL 3.2.2 from homebrew)

#### Supported Satellites and Sensors

* Sentinel-1A/B

## Installation

An installation guide can be found at [insar/docs/Installation.md](insar/docs/Installation.md)

## GAMMA-INSAR Unit Testing

Running unit tests for `gamma_insar` is as simple as running `pytest` from the project directory on a supported platform (docker options below).  The test suite was written with the assumption that GAMMA is *unavailable*, but all other dependencies are required - this is suitable for testing on a wider range of systems (such as developers machines & CI testing environments which likely won't have GAMMA licenses).

As a result of this design decision, the *unit* tests only test the logic of the workflow - not the correctness of the processed data.

To run unit tests:
```BASH
cd ~/gamma_insar
source configs/activateNCI.env ~/gamma_insar_install
pytest --disable-warnings -q  # should be error free
```

Code coverage can be checked with pytest-cov. Note that running `coverage.py` alone with this repo **does not accurately record coverage results!** The `pytest-cov` tool is required to measure coverage correctly.

To measure code test coverage:

```BASH
# run tests & display coverage report at the terminal
pytest -q --disable-warnings --cov=insar

# run tests & generate an interactive HTML report
pytest -q --disable-warnings --cov-report=html --cov=insar
```

The report is saved to `coverage_html_report` in the project dir.

The unit tests may also be run on platforms which are unsupported, or do not have the dependencies installed - via docker (assuming the target platform can run docker images), please refer to the comments at the top of the `Dockerfile` for instructions on building the image and running `pytest` in that image.

### The InSAR Workflow

Several preliminary steps are required before the main processing workflow can be executed, these steps produce a database of data acquisitions that 'can' be processed - which is queried by the workflow to determine what data to process for a proided ESRI shapefile (to geometrically bound the area of interest) and date-range.

These steps are:
1) [Extract metadata](#Metadata-YAML-extraction) from the source data (eg: Sentinel-1 SLC zip files) and store
   that information into yaml files. This needs to be done for a time
   period any time data has been added/removed/changed.

2) [Generate a sqlite database](#Database-Creation) from all the yamls created in (1).
   This db file is used by the workflow to quickly query/filter scene information for processing.

3) [Produce a set of shapefiles](#shapefile-creation) that cover the areas of interest for processing.

4) [Process data products](#Data-processing) (eg: backscatter and/or interferograms)

5) [Package data products](#Product-packaging) for ODC indexing

### Documentation

Various documents are available to help users and developers learn more about the project, what it does, how to use it, and it's technical structure.

User guides:
 * [Installation guide](insar/docs/Installation.md)
 * [Running the workflow](insar/docs/RunningTheWorkflow.md)
 * [Geospatial database indexing](insar/docs/DatabaseIndexing.md)
 * [Data packaging](insar/docs/Packaging.md)
 * [Example .proc file](template.proc)

Technical Documents:
 * [InSAR Stack description](insar/docs/Stack.md)
 * [Metadata documentation](insar/docs/Metadata.md)

Developer guides:
 * [InSAR workflow overview](insar/docs/Workflow.md)

### Unit testing on NCI

The Gadi system at NCI is a supported platform, and as such running unit tests is no different - simply enter/activate your installation (see [Installation on NCI](#Installation-on-NCI)) via `configs/activateNCI.env` and run `pytest` from the project directory.


### Shapefile creation

Shapefiles should typically match data acquisition patterns for the sensor data being processed, this is because mismatches in available data across acquisitions can cause complex problems in `gamma_insar` due to requirements imposed by the coregistration & the interferogram tree linking algorithm.

Due to how coregistration & interferogram trees work - only dates that share identical "sets geographic regions" in the source data (eg: bursts in Sentinel-1) may be correlated, and thus any acquisition that does **not** share the most expansive set of data (eg: is missing a burst) will be excluded.

The reason for this largely comes down to the fact that both coregistration and interferograms are in their nature an operation between two distinct scenes, and thus if data in scene A does not exist in scene B there is nothing to coregister with nor produce an interferogram from...

For this reason it is strongly recommended that shapefiles are produced in such a way that all scenes will consistently have the same set of data.
