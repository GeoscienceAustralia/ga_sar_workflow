#!/bin/bash

display_usage() {
    echo ""
    echo "*******************************************************************************"
    echo "* coregister_S1_slave_SLC: Coregisters Sentinel-1 IWS SLC to chosen master    *"
    echo "*                          SLC geometry                                       *"
    echo "*                                                                             *"
    echo "* input:  [proc_file]  name of GAMMA proc file (eg. gamma.proc)               *"
    echo "*         [slave]      slave scene ID (eg. 20120520)                          *"
    echo "*         [rlks]       range multi-look value (for SLCs: from *.mli.par file  *"
    echo "*                      or for ifms: from proc file)                           *"
    echo "*         [alks]       azimuth multi-look value (for SLCs: from *.mli.par     *"
    echo "*                      file or for ifms: from proc file)                      *"
    echo "*                                                                             *"
    echo "* author: Matt Garthwaite @ GA       12/05/2015, v1.0                         *"
    echo "*******************************************************************************"
    echo -e "Usage: coregister_S1_slave_SLC.bash [proc_file] [slave] [rlks] [alks]"
    }

if [ $# -lt 4 ]
then 
    display_usage
    exit 1
fi

proc_file=$1
slave=$2
rlks=$3
alks=$4

## Variables from parameter file (*.proc)
platform=`grep Platform= $proc_file | cut -d "=" -f 2`
project=`grep Project= $proc_file | cut -d "=" -f 2`
sensor=`grep Sensor= $proc_file | cut -d "=" -f 2`
track_dir=`grep Track= $proc_file | cut -d "=" -f 2`
master=`grep Master_scene= $proc_file | cut -d "=" -f 2`
polar=`grep Polarisation= $proc_file | cut -d "=" -f 2`
subset=`grep Subsetting= $proc_file | cut -d "=" -f 2`
subset_done=`grep Subsetting_done= $proc_file | cut -d "=" -f 2`
snr=`grep coreg_snr_thresh= $proc_file | cut -d "=" -f 2`
npoly=`grep coreg_model_params= $proc_file | cut -d "=" -f 2`
win=`grep coreg_window_size= $proc_file | cut -d "=" -f 2`
nwin=`grep coreg_num_windows= $proc_file | cut -d "=" -f 2`
ovr=`grep coreg_oversampling= $proc_file | cut -d "=" -f 2`
niter=`grep coreg_num_iterations= $proc_file | cut -d "=" -f 2`

## Identify project directory based on platform
if [ $platform == NCI ]; then
    proj_dir=/g/data1/dg9/INSAR_ANALYSIS/$project/$sensor/GAMMA
else
    proj_dir=/nas/gemd/insar/INSAR_ANALYSIS/$project/$sensor/GAMMA
fi

slc_dir=$proj_dir/$track_dir/`grep SLC_dir= $proc_file | cut -d "=" -f 2`
dem_dir=$proj_dir/$track_dir/`grep DEM_dir= $proc_file | cut -d "=" -f 2`

cd $proj_dir

## Insert scene details top of NCI .e file
echo "" 1>&2 # adds spaces at top so scene details are clear
echo "" 1>&2
echo "PROCESSING_SCENE: "$project $track_dir $slave $rlks"rlks" $alks"alks" 1>&2

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

master_dir=$slc_dir/$master
slave_dir=$slc_dir/$slave
master_slc_name=$master"_"$polar
slave_slc_name=$slave"_"$polar
master_mli_name=$master"_"$polar"_"$rlks"rlks"
slave_mli_name=$slave"_"$polar"_"$rlks"rlks"

## files located in SLC directories
master_mli=$master_dir/r$master_mli_name.mli
master_mli_par=$master_mli.par
master_slc=$master_dir/r$master_slc_name.slc
master_slc_par=$master_slc.par
master_slc_tab=$master_dir/slc_tab

slave_mli=$slave_dir/$slave_mli_name.mli
slave_mli_par=$slave_mli.par
slave_slc=$slave_dir/$slave_slc_name.slc
slave_slc_par=$slave_slc.par
slave_slc_tab=$slave_dir/slc_tab

slc1=r$slave_slc_name"_IW1.slc"
slc1_par=$slc1.par
tops_par1=$slc1.TOPS_par
slc2=r$slave_slc_name"_IW2.slc"
slc2_par=$slc2.par
tops_par2=$slc2.TOPS_par
slc3=r$slave_slc_name"_IW3.slc"
slc3_par=$slc3.par
tops_par3=$slc3.TOPS_par

rslc=$slave_dir/r$slave_slc_name.slc 
rslc_par=$rslc.par
rmli=$slave_dir/r$slave_mli_name.mli 
rmli_par=$rmli.par
rslc_tab=$slave_dir/rslc_tab

lt=$slave_dir/$master-$slave_mli_name.lt
off=$slave_dir/$master-$slave_mli_name.off

## Set up coregistration results file
check_file=$slave_dir/slave_coregistration_results"_"$rlks"_rlks_"$alks"_alks.txt"
#if [ -f $check_file ]; then
#    rm -f $check_file 
#else
#    :
#fi
echo "Slave_Coregistration_Results_"$rlks"_rlks_"$alks"_alks" > $check_file
echo "final model fit std. dev. (samples)" >> $check_file
echo "Ref Master" > temp1_$rlks
echo "Slave" > temp2_$rlks
echo "Range" > temp3_$rlks
echo "Azimuth" > temp4_$rlks
paste temp1_$rlks temp2_$rlks temp3_$rlks temp4_$rlks >> $check_file
rm -f temp1_$rlks temp2_$rlks temp3_$rlks temp4_$rlks

cd $slave_dir

## Determine range and azimuth looks in MLI
echo " "
echo "MLI range and azimuth looks: "$rlks $alks
echo " "

#-------------------------

## files located in DEM directory
rdc_dem=$dem_dir/$master_mli_name"_rdc.dem"

## Generate initial lookup table between master and slave MLI considering terrain heights from DEM coregistered to master
GM rdc_trans $master_mli_par $rdc_dem $slave_mli_par lt0

slave_mli_width=`awk 'NR==11 {print $2}' $slave_mli_par`
master_mli_width=`awk 'NR==11 {print $2}' $master_mli_par`
slave_mli_length=`awk 'NR==12 {print $2}' $slave_mli_par`

GM geocode lt0 $master_mli $master_mli_width $rmli $slave_mli_width $slave_mli_length 2 0

GM create_diff_par $slave_mli_par $slave_mli_par diff.par 1 0

## Measure offset between slave MLI and resampled slave MLI
GM init_offsetm $rmli $slave_mli diff.par 1 1

GM offset_pwrm $rmli $slave_mli diff.par offs0 snr0 - - - 2

## Fit the offset only
GM offset_fitm offs0 snr0 diff.par coffs0 - 7.0 1

## Refinement of initial geocoding look up table
GM gc_map_fine lt0 $master_mli_width diff.par $lt

## Create table for resampled burst SLCs
rm -f $rslc_tab
for swath in 1 2 3; do
    bslc="slc$swath"
    bslc_par=${!bslc}.par
    btops="tops_par$swath"
    echo $slave_dir/${!bslc} $slave_dir/$bslc_par $slave_dir/${!btops} >> $rslc_tab
done

## Resample slave SLC into geometry of master SLC using lookup table and generate mosaic SLC    
GM SLC_interp_lt_S1_TOPS $slave_slc_tab $slave_slc_par $master_slc_tab $master_slc_par $lt $master_mli_par $slave_mli_par - $rslc_tab $rslc $rslc_par

#------------------------

## set up iterable loop
i=1
while [ $i -le $niter ]; do

    ioff=$off$i
    rm -f offs snr offsets coffsets
    echo "Starting Iteration "$i

## Measure offsets for refinement of lookup table using initially resampled slave SLC
    GM create_offset $master_slc_par $rslc_par $ioff 1 $rlks $alks 0

## No SLC oversampling for S1 due to strong Doppler centroid variation in azimuth
    GM offset_pwr $master_slc $rslc $master_slc_par $rslc_par $ioff offs snr 256 64 offsets $ovr $nwin $nwin $snr 

## Fit constant offset term only for S1 due to short length orbital baselines
    GM offset_fit offs snr $ioff - coffsets 10.0 $npoly 0

## Create blank offset file for first iteration and calculate the total estimated offset
    if [ $i == 1 ]; then
	GM create_offset $master_slc_par $rslc_par $off"0" 1 $rlks $alks 0

	GM offset_add $off"0" $ioff $off
    else
## Calculate the cumulative total estimated offset
	GM offset_add $off $ioff $off
    fi

## if azimuth offset is less than 0.02 and range offset is less than 0.2 then break iterable loop. Precision azimuth coregistration is essential for S1 IWS mode interferometry
    azoff=`grep "final azimuth offset poly. coeff." output.log | tail -2 | head -1 | awk '{print $6}'`
    rgoff=`grep "final range offset poly. coeff." output.log | tail -2 | head -1 | awk '{print $6}'`  
    test1=`echo $azoff | awk '{if ($1 < 0) $1 = -$1; printf "%i\n", $1*100}'`
    test2=`echo $rgoff | awk '{if ($1 < 0) $1 = -$1; printf "%i\n", $1*10}'`
    echo "Iteration "$i": azimuth offset is "$azoff", range offset is "$rgoff

## Perform resampling of slave SLC using lookup table and offset model, and generate mosaic SLC
    GM SLC_interp_lt_S1_TOPS $slave_slc_tab $slave_slc_par $master_slc_tab $master_slc_par $lt $master_mli_par $slave_mli_par $off $rslc_tab $rslc $rslc_par

    if [ $test1 -lt 2 -a $test2 -lt 2 ]; then
	break
    fi
    i=$(($i+1))
done

#-------------------------

## Determine a refinement to the azimuth offset estimation in the burst overlap regions
GM S1_coreg_overlap $master_slc_tab $rslc_tab $master-$slave $off $off".corrected" 0.8 100

## Perform fourth resampling of slave SLC using lookup table and corrected offset information
GM SLC_interp_lt_S1_TOPS $slave_slc_tab $slave_slc_par $master_slc_tab $master_slc_par $lt $master_mli_par $slave_mli_par $off".corrected" $rslc_tab $rslc $rslc_par

#-------------------------

GM multi_look $rslc $rslc_par $rmli $rmli_par $rlks $alks

rm -f offs0 snr0 coffs0 offs snr coffs coffsets lt0

## Extract final model fit values to check coregistration
echo $master > temp1_$rlks
echo $slave > temp2_$rlks
grep "final" temp2 > temp3_$rlks
awk '{print $8}' temp3_$rlks > temp4_$rlks
awk '{print $10}' temp3_$rlks > temp5_$rlks
paste temp1_$rlks temp2_$rlks temp4_$rlks temp5_$rlks >> $check_file
rm -f temp*


# script end 
####################

## Copy errors to NCI error file (.e file)
if [ $platform == NCI ]; then
   cat error.log 1>&2
#   rm temp_log
else
    $slave_dir/temp_log
fi
