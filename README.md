## GAMMA-INSAR

A tool to process Sentinel-1 SLC to Aanalysis Ready Data using GAMMA SOFTWARE.

## Installation

    python setup.py install --prefix=<prefix> 

Python 3.6+ is supported.

## Operating System tested
Linux

## Supported Satellites and Sensors
* Sentinel-1A/B

## Requirements
* [attrs>=17.4.0]
* [Click>=7.0]
* [GDAL>=2.4]
* [geopandas>=0.4.1]
* [luigi>=2.8.3]
* [matplotlib>=3.0.3]
* [numpy>=1.8]
* [pandas>-0.24.2]
* [pyyaml>=3.11]
* [rasterio>=1,!=1.0.3.post1,!=1.0.3]
* [structlog>=16.1.0]
* [shapely>=1.5.13]
* [spatialist==0.4]
* [GAMMA-SOFTWARE >= June 2019 release]

`export PYTHONPATH=<path-to-gamma-software>:$PYTHONPATH`

## Usage

`gamma_insar': Process SLC data to ARD from the commandline.

	$gamma_insar ARD --help

	usage: gamma_insar ARD
		   [REQUIRED PARAMETERS]
		   --vector-file VECTOR_FILE		A full path to a Sentinel-1 tract and frame vector-file 
		   --start-date START_DATE		A start-date['YYYY-MM-DD'] of SLC data acquisition
		   --end-date END_DATE			An end-date['YYYY-MM-DD'] of SLC data acquisition
		   --workdir WORKDIR			A full path to a working directory to output logs
	           --outdir OUTDIR	 		A full path to an output directory
		   --polarization POLARIZATION  	Polarizations to be processed [VV|VH]	
		   --cleanup CLEANUP			A Flag[yes|no] to specify to clean up intermediary files 
							Highly recommended to cleanup
		   --database-name DATABASE_NAME	A full path to SLC-metata database with burst informations
		   --orbit ORBIT			A Sentinel-1 orbit [A|D]
		   --dem-img DEM_IMG			A full path to a Digital Elevation Model
		   --multi-look MULTI_LOOK		A multi-look value
		   --poeorb-path POEORB_PATH		A full path to a directory with precise orbit file
		   --resorb-path RESORB_PATH		A full path to a directory with restitution orbit file
		   --num-threads NUM_THREADS		A number of threads to be used Environmental variable 
							OMP_NUM_THREADS gets modified while GAMMA SOFTWARE is called 
							for co-registration and inteferograms processing.
		   --workers WORKERS			Number of workers assigned to a luigi scheduler 
		   --local-scheduler SCHEDULER		Use only local-scheduler


`gamma_insar ARD --vector-file <path-to-vector-file> --start-date <start-date> --end-date <end-date> --workdir <path-to-workdir> --outdir <path-to-outdir> --workers <number-of-workers> --local-scheduler` 




