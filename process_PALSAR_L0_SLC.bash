#!/bin/bash

display_usage() {
    echo ""
    echo "*******************************************************************************"
    echo "* process_PALSAR_L0_SLC:  Script takes Level 1.0 (raw) PALSAR1 fine beam data *"
    echo "*                         and produces sigma0 calibrated SLC and and Multi-   *"
    echo "*                         look intensity (MLI) images in GAMMA format.        *"
    echo "*                                                                             *"
    echo "*                         Requires a 'frame.list' text file to be created in  *"
    echo "*                         the project directory. This lists the frame numbers *"
    echo "*                         on each line (e.g. 7160).                           *"
    echo "*                                                                             *"
    echo "* input:  [proc_file]    name of GAMMA proc file (eg. gamma.proc)             *"
    echo "*         [scene]        scene ID (eg. 20070112)                              *"
    echo "*         [rlks]         MLI range looks                                      *"
    echo "*         [alks]         MLI azimuth looks                                    *"
    echo "*                                                                             *"
    echo "* author: Sarah Lawrie @ GA       01/05/2015, v1.0                            *"
    echo "*******************************************************************************"
    echo -e "Usage: process_PALSAR_L0_SLC.bash [proc_file] [scene] [rlks] [alks]"
    }

if [ $# -lt 4 ]
then 
    display_usage
    exit 1
fi

if [ $2 -lt "10000000" ]; then 
    echo "ERROR: Scene ID needed in YYYYMMDD format"
    exit 1
else
    scene=$2
fi

proc_file=$1

## Variables from parameter file (*.proc)
platform=`grep Platform= $proc_file | cut -d "=" -f 2`
project=`grep Project= $proc_file | cut -d "=" -f 2`
track_dir=`grep Track= $proc_file | cut -d "=" -f 2`
polar=`grep Polarisation= $proc_file | cut -d "=" -f 2`
subset=`grep Subsetting= $proc_file | cut -d "=" -f 2`
sensor=`grep Sensor= $proc_file | cut -d "=" -f 2`
frame_list=`grep List_of_frames= $proc_file | cut -d "=" -f 2`
raw_dir_ga=`grep Raw_data_GA= $proc_file | cut -d "=" -f 2`
raw_dir_mdss=`grep Raw_data_MDSS= $proc_file | cut -d "=" -f 2`

slc_rlks=$3
slc_alks=$4


## Identify project directory based on platform
if [ $platform == NCI ]; then
    proj_dir=/g/data1/dg9/INSAR_ANALYSIS/$project/$sensor/GAMMA
    raw_dir=$proj_dir/raw_data/$track_dir
else
    proj_dir=/nas/gemd/insar/INSAR_ANALYSIS/$project/$sensor/GAMMA
    raw_dir=$raw_dir_ga
fi

cd $proj_dir/$track_dir

## Insert scene details top of NCI .e file
echo "" 1>&2 # adds spaces at top so scene details are clear
echo "" 1>&2
echo "PROCESSING_PROJECT: "$project $track_dir $scene 1>&2

## Copy output of Gamma programs to log files
GM()
{
    echo $* | tee -a command.log
    echo
    $* >> output.log 2> temp_log
    cat temp_log >> error.log
    #cat output.log (option to add output results to NCI .o file if required)
}

## Load GAMMA based on platform
if [ $platform == NCI ]; then
    GAMMA=`grep GAMMA_NCI= $proc_file | cut -d "=" -f 2`
    source $GAMMA
else
    GAMMA=`grep GAMMA_GA= $proc_file | cut -d "=" -f 2`
    source $GAMMA
fi

slc_dir=$proj_dir/$track_dir/`grep SLC_dir= $proc_file | cut -d "=" -f 2`
scene_dir=$slc_dir/$scene

echo " "
echo "MLI range and azimuth looks: "$slc_rlks $slc_alks
echo " "

mkdir -p $slc_dir
cd $slc_dir
mkdir -p $scene
cd $scene_dir

raw_file_list=raw_file_list
rm -f $scene_dir/$raw_file_list 

## File names
slc_name=$scene"_"$polar
mli_name=$scene"_"$polar"_"$slc_rlks"rlks"
sensor_antpat=$MSP_HOME/sensors/palsar_ant_20061024.dat
msp_antpat="PALSAR_antpat_MSP_"$polar.dat
sensor_par="PALSAR_sensor_"$polar.par
msp_par=p$slc_name.slc.par
para=$slc_name"_SLC_parameters.txt"
raw=$slc_name.raw
slc=$slc_name.slc
slc_par=$slc.par
mli=$mli_name.mli
mli_par=$mli.par
tiff=$mli_name.tif
ras_out=$mli_name.ras
fbd2fbs_slc=$slc_name"_FBS.slc"
fbd2fbs_par=p$slc_name"_FBS.slc.par"

## Set mode based on polarisation
pol_list=$scene_dir/pol_list
rm -f $pol_list

if [ ! -e $slc_dir/$scene/$slc ]; then
    while read frame_num; do
	if [ ! -z $frame_num ]; then
	    frame=`echo $frame_num | awk '{print $1}'`
	    ls $raw_dir/F$frame/date_dirs/$scene/IMG-HH* >& hh_temp
	    ls $raw_dir/F$frame/date_dirs/$scene/IMG-HV* >& hv_temp
	    temp="ls: cannot access"
	    temp1=`awk '{print $1" "$2" "$3}' hh_temp`
	    if [ "$temp1" == "$temp" ]; then
		:
	    else
		basename $raw_dir/F$frame/date_dirs/$scene/IMG-HH* >> $pol_list 
	    fi
	    temp2=`awk '{print $1" "$2" "$3}' hv_temp`
	    if [ "$temp2"  == "$temp" ]; then
		:
	    else
		basename $raw_dir/F$frame/date_dirs/$scene/IMG-HV* >> $pol_list
	    fi
	fi
	rm -rf hh_temp hv_temp
    done < $proj_dir/$track_dir/$frame_list

    num_hv=`grep -co "HV" $pol_list`
    if [ "$num_hv" -eq 0 -a "$polar" == HH ]; then 
	mode=FBS
    elif [ "$num_hv" -ge 1 -a "$polar" == HH ]; then 
	mode=FBD
    elif [ $polar == HV ]; then
	mode=FBD
    else
	:
    fi
    echo "Mode:" $mode "  Polarisation:" $polar
    rm -f $pol_list

## Produce raw data files
    while read frame_num; do
	if [ ! -z $frame_num ]; then
	    frame=`echo $frame_num | awk '{print $1}'`
	    ls $raw_dir/F$frame/date_dirs/$scene/LED-ALP* >& temp
	    LED=`awk '{print $1}' temp`
	    rm -f temp
            # Check polarisation
	    if [ $polar == HH ]; then
		ls $raw_dir/F$frame/date_dirs/$scene/IMG-HH* >& temp
		IMG=`awk '{print $1}' temp`
		rm -f temp
	    else 
		ls $raw_dir/F$frame/date_dirs/$scene/IMG-HV* >& temp
		IMG=`awk '{print $1}' temp`
		rm -f temp
	    fi
	    sensor_fm_par="PALSAR_sensor_"$frame"_"$polar.par
	    msp_fm_par="p"$scene"_"$frame"_"$polar.slc.par
	    raw_fm_file=$scene"_"$frame"_"$polar.raw

            # Set polarisation (0 = H, 1 = V)
	    tx_pol=`echo $polar | awk '{print substr($1,1,1)}'`
	    rx_pol=`echo $polar | awk '{print substr($1,2,1)}'`
            if [ $tx_pol == H ]; then
		tx_pol1=0
	    else
		tx_pol1=1
	    fi
	    if [ $rx_pol == H ]; then
		rx_pol1=0
	    else
		rx_pol1=1
	    fi
            ## Generate the processing parameter files and raw data in GAMMA format
	    GM PALSAR_proc $LED $sensor_fm_par $msp_fm_par $IMG $raw_fm_file $tx_pol1 $rx_pol1
            ## Copy raw data file details to text file to check if concatenation of scenes along track is required
	    echo $raw_fm_file $sensor_fm_par $msp_fm_par >> $raw_file_list
	fi
    done < $proj_dir/$track_dir/$frame_list

    ## Check if scene concatenation is required (i.e. a scene has more than one frame)
    lines=`awk 'END{print NR}' $raw_file_list`
    if [ $lines -eq 1 ]; then
        ## rename files to enable further processing (remove reference to 'frame' in file names)
	mv $sensor_fm_par $sensor_par
	mv $msp_fm_par $msp_par
	mv $raw_fm_file $raw
	rm -f $raw_file_list
    else
        ## Concatenate scenes into one output raw data file
	GM cat_raw $raw_file_list $sensor_par $msp_par $raw 1 0 - 
	rm -f p$scene"_"*"_"$polar.slc.par
	rm -f "PALSAR_sensor_"*"_"$polar.par
	rm -f $scene"_"*"_"$polar.raw
	rm -f $raw_file_list
    fi

    ## Correction for antenna pattern
    GM PALSAR_antpat $sensor_par $msp_par $sensor_antpat $msp_antpat - $tx_pol1 $rx_pol1

    ## Determine the Doppler Ambiguity
    ## Use dop_mlcc instead of dop_ambig when number of raw echoes greater than 8192
    GM dop_mlcc $sensor_par $msp_par $raw $slc_name.mlcc
    plot_mlcc.bash $slc_name.mlcc

    ## Estimate the doppler centroid with cross correlation method
    ## If result of 'doppler' shows that linear model is not good enough use 'azsp_IQ' to determine constant value
    GM doppler $sensor_par $msp_par $raw $slc_name.dop
    plot_dop.bash $slc_name.dop

    ## Estimate the range power spectrum
    ## Look for potential radio frequency interference (RFI) to the SAR signal
    GM rspec_IQ $sensor_par $msp_par $raw $slc_name.rspec

    ## Check range spectrum for spikes indicating RFI. If they exist can be suppresed during range compression 'pre_rc'
    #plot_rspec.bash $slc_name.rspec $scene $raw_dir $msp_par

    ## Range compression
    ## second to last parameter is for RFI suppression.
    GM pre_rc $sensor_par $msp_par $raw $slc_name.rc - - - - - - - - 0 -

    ## REMOVE RAW IMAGE FILE HERE
    rm -f $raw

    ## Autofocus estimation and Azimuth compression (af replaces autof)
    ## run az_proc and af twice, DO NOT run af mutliple times before reprocessing image
    ## default SNR threshold is 10
    ## slc calibrated as sigma0 (assume 34.3 angle)
    if [ $polar == HH -a $mode == FBS ]; then
	cal_const=-51.9
    elif [ $polar == HH -a $mode == FBD ]; then
	cal_const=-51.8
    elif [ $polar == HV -a $mode == FBD ]; then
	cal_const=-58.3
    else
	:
    fi
    GM az_proc $sensor_par $msp_par $slc_name.rc $slc 8192 0 $cal_const 0
    GM af $sensor_par $msp_par $slc 1024 4096 - - 10 1 0 0 $slc_name.af
    GM az_proc $sensor_par $msp_par $slc_name.rc $slc 8192 0 $cal_const 0
    GM af $sensor_par $msp_par $slc 1024 4096 - - 10 1 0 0 $slc_name.af

    ## REMOVE RC FILE HERE
    rm -f $slc_name.rc

    ## Generate ISP SLC parameter file from MSP SLC parameter file (for full SLC)
    GM par_MSP $sensor_par $msp_par $slc_par 0

    ## FBD to FBS Conversion
    if [ $polar == HH -a $mode == FBD ]; then
	GM SLC_ovr $slc $slc_par $fbd2fbs_slc $fbd2fbs_par 2
	rm -f $slc
	rm -f $slc_par
	mv $fbd2fbs_slc $slc
	mv $fbd2fbs_par $slc_par
    else
	:
    fi

## Add details from ISP parameter file to MSP parameter file
#temp1=`awk 'NR==4 {print $2}' $slc_par`
#sed -i "4i sensor:                 $temp1" $msp_par 
#temp2=`grep Track= $proc_file | cut -d "=" -f 2`
#sed -i "6i track:                  $temp2" $msp_par 

              #### sort frame information for later - need to incorporate concatenation
#temp3=`echo $raw_dir | cut -d "_" -f7`
#temp3=9999
#sed -i "7i frame:                  $temp3" $msp_par 
              #### no orbit information
#temp4=`grep Orbit: $raw_report | awk '{print $4}'`
#temp4=-
#sed -i "8i orbit:                  $temp4" $msp_par 
#temp5=`grep Orientation= $proc_file | cut -d "=" -f 2`
#sed -i "9i orientation:            $temp5" $msp_par 
#temp6=`ls -ltrh $slc | awk '{print $5}'`
#sed -i "10i slc_size:               $temp6" $msp_par 
#temp7=`grep start_time: $slc_par | awk '{print $2}'`
#sed -i "12i start_time:             $temp7" $msp_par 
#temp8=`grep center_time: $slc_par | awk '{print $2}'`
#sed -i "13i center_time:            $temp8" $msp_par 
#temp9=`grep end_time: $slc_par | awk '{print $2}'`
#sed -i "14i end_time:               $temp9" $msp_par 
#temp10=`grep incidence_angle: $slc_par | awk '{print $2}'`
#sed -i "16i incidence_angle:        $temp10" $msp_par 
#temp11=`grep radar_frequency: $slc_par | awk '{print $2}'`
#sed -i "17i radar_frequency:        $temp11" $msp_par 
#temp12=L
#sed -i "18i sar_band:               $temp12" $msp_par 
#temp13=`grep sar_to_earth_center: $slc_par | awk '{print $2}'`
#sed -i "21i sar_to_earth_center:          $temp13" $msp_par 
#temp14=`grep earth_radius_below_sensor: $slc_par | awk '{print $2}'`
#sed -i "22i earth_radius_below_sensor:    $temp14" $msp_par 
#temp15=`grep first_slant_range_polynomial: $slc_par | awk '{print $2}'`
#sed -i "42i first_slant_range_polynomial:      $temp15" $msp_par 
#temp16=`grep center_slant_range_polynomial: $slc_par | awk '{print $2}'`
#sed -i "43i center_slant_range_polynomial:          $temp16" $msp_par 
#temp17=`grep last_slant_range_polynomial: $slc_par | awk '{print $2}'`
#sed -i "44i last_slant_range_polynomial:   $temp17" $msp_par 
#temp18=`grep azimuth_line_time: $slc_par | awk '{print $2}'`
#sed -i "75i azimuth_line_time:                $temp18" $msp_par 
#temp19=`grep azimuth_angle: $slc_par | awk '{print $2}'`
#sed -i "76i azimuth_angle:                  $temp19" $msp_par 
#temp20=`grep azimuth_scale_factor: $slc_par | awk '{print $2}'`
#sed -i "77i azimuth_scale_factor:             $temp20" $msp_par 
#temp21=`grep range_scale_factor: $slc_par | awk '{print $2}'`
#sed -i "78i range_scale_factor:               $temp21" $msp_par 
#temp22=`grep line_header_size: $slc_par | awk '{print $2}'`
#sed -i "81i line_header_size:                 $temp22" $msp_par 
#temp23=`grep adc_sampling_rate: $slc_par | awk '{print $2}'`
#sed -i "82i adc_sampling_rate:                 $temp23" $msp_par 
#temp24=`grep chirp_bandwidth: $slc_par | awk '{print $2}'`
#sed -i "83i chirp_bandwidth:                  $temp24" $msp_par 
#temp25=`grep image_geometry: $slc_par | awk '{print $2}'`
#sed -i "85i image_geometry:                   $temp25" $msp_par 

## Create copy of MSP parameter file to create generic SLC parameter file
#cp $msp_par $para
#sed -i "1s/.*/SLC_Parameter_File/" $para

else
    echo " "
    echo "Full SLC already created."
    echo " "
fi

## Multi-look SLC
GM multi_look $slc $slc_par $mli $mli_par $slc_rlks $slc_alks 0


## Create low-res preview tiff
#mli_width=`grep range_samples: $mli_par | awk '{print $2}'`
#GM data2tiff $mli $mli_width 2 $tiff

## Create low-res ras image (for location plot)
#GM raspwr $mli $mli_width 1 0 1 1 1 0.35 1 $ras_out 0 0

## corner coordinates given in SLC MSP parameter file
#grep map_coordinate_4 $msp_par | awk '{print $2, $3}' > slc_coords
#grep map_coordinate_2 $msp_par | awk '{print $2, $3}' >> slc_coords
#grep map_coordinate_1 $msp_par | awk '{print $2, $3}' >> slc_coords
#grep map_coordinate_3 $msp_par | awk '{print $2, $3}' >> slc_coords
#grep map_coordinate_5 $msp_par | awk '{print $2, $3}' >> slc_coords

## Make SLC location plot
#plot_SLC_loc.bash $proc_file $scene $msp_par $sensor $ras_out



# script end 
####################

## Copy errors to NCI error file (.e file)
if [ $platform == NCI ]; then
   cat error.log 1>&2
   rm temp_log
else
   rm temp_log
fi
