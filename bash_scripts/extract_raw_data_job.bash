#!/bin/bash

display_usage() {
    echo ""
    echo "*******************************************************************************"
    echo "* extract_raw_data_job: scripts used to extract raw data                      *"
    echo "*                                                                             *"
    echo "*                                                                             *"
    echo "* input:  [proc_file]   name of GAMMA proc file (eg. gamma.proc)              *"
    echo "*                                                                             *"
    echo "* author: Sarah Lawrie @ GA       13/08/2018, v1.0                            *"
    echo "*             							                                    *"
    echo "*******************************************************************************"
    echo -e "Usage: extract_raw_data_job.bash [proc_file]"
    }

if [ $# -lt 1 ]
then
    display_usage
    exit 1
fi 
proc_file=$1

##########################   GENERIC SETUP  ##########################

# Load generic GAMMA functions
source ~/repo/gamma_insar/gamma_functions

# Load variables and directory paths
proc_variables $proc_file
final_file_loc
# Load GAMMA to access GAMMA programs
source $config_file

# Print processing summary to .o & .e files
PBS_processing_details $project $track $scene

######################################################################

pbs_job_dirs
final_file_loc
proc_file="$(basename -- $proc_file)"

if [ $do_raw == yes ]; then
    cd $extract_raw_batch_dir
    if [ -e all_raw_job_ids ]; then
	    echo "Raw data already extracted."
	    echo ""
    else
	    echo ""
	    echo "Extracting raw data ..."
	    rm -f list # temp file used to collate all PBS job numbers to dependency list

        if [ $sensor == 'S1' ]; then
            queue1=$queue
            input_list=$s1_file_list
            awk '/FILES_TO_DOWNLOAD/ { show=1 } show; /SUBSET_BURSTS/ { show=0 }' $input_list | tail -n+3 | head -n -2 > $list_dir/download_list
            list=$list_dir/download_list
            nlines=`cat $list | sed '/^\s*$/d' | wc -l`

            # create frame subset list (cut scenes to frame extents)
            awk '/SUBSET_BURSTS/ { show=1 } show; /ORG_BURSTS_V_MASTER_BURSTS/ { show=0 }' $input_list | tail -n+3 | head -n -2 | awk '{print $1,$5,$6,$7,$8}' > $list_dir/temp
            if [ -e $list_dir/frame_subset_list ]; then
                rm -rf $list_dir/frame_subset_list
            fi
            while read subset; do
                date=`echo $subset | awk '{print $1}'`
                start_iw1=`echo $subset | awk '{print $2}' | cut -d '-' -f 1`
                stop_iw1=`echo $subset | awk '{print $2}' | cut -d '-' -f 2`
                start_iw2=`echo $subset | awk '{print $3}' | cut -d '-' -f 1`
                stop_iw2=`echo $subset | awk '{print $3}' | cut -d '-' -f 2`
                start_iw3=`echo $subset | awk '{print $4}' | cut -d '-' -f 1`
                stop_iw3=`echo $subset | awk '{print $4}' | cut -d '-' -f 2`
                complete_frame=`echo $subset | awk '{print $5}' | cut -d '-' -f 1`
                echo $date "1" $start_iw1 $stop_iw1 $complete_frame >> $list_dir/frame_subset_list
                echo $date "2" $start_iw2 $stop_iw2 $complete_frame >> $list_dir/frame_subset_list
                echo $date "3" $start_iw3 $stop_iw3 $complete_frame >> $list_dir/frame_subset_list
                done < $list_dir/temp
                rm -rf $list_dir/temp
        else
            queue1=$mdss_queue
            list=$scene_list
            nlines=`cat $list | sed '/^\s*$/d' | wc -l`
        fi

        echo Need to process $nlines files
        # PBS parameters
        wt1=`echo $raw_walltime | awk -F: '{print ($1*60) + $2 + ($3/60)}'` # walltime for a single process_slc in minutes
        pbs_job_prefix=raw_
        script=extract_raw_data.bash
        script_type=-
        depend_job=0 #no dependencies
        depend_type=-
        job_type=1 #1 for batch job, 2 for manual job

        # Work out number of jobs to run within maximum number of jobs allowed and create jobs
        if [ $nlines -le $minjobs ]; then
            jobs1=$nlines
            steps1=1
            jobs2=0
            steps2=0
        else
            steps2=$((nlines/minjobs))
            steps1=$((nlines%minjobs))
            jobs1=$steps1
            steps1=$((steps2+1))
            jobs2=$((minjobs-jobs1))
        fi

        echo Preparing to run $jobs1 jobs with $steps1 steps and $jobs2 jobs with $steps2 steps processing $((jobs1*steps1+jobs2*steps2)) files
        j=0
        {
            multi_jobs $pbs_run_loc $pbs_job_prefix $nci_project $raw_mem $raw_ncpus $queue1 $script $depend_job $depend_type $job_type $extract_raw_batch_dir $script_type jobs1 steps1 j
            multi_jobs $pbs_run_loc $pbs_job_prefix $nci_project $raw_mem $raw_ncpus $queue1 $script $depend_job $depend_type $job_type $extract_raw_batch_dir $script_type jobs2 steps2 jobs1
        } < $list

        # Create manual PBS jobs
        cd $extract_raw_manual_dir
        job_type=2 #1 for batch job, 2 for manual job
        depend_job=0
        j=0
        {
            multi_jobs $pbs_run_loc $pbs_job_prefix $nci_project $raw_mem $raw_ncpus $queue1 $script $depend_job $depend_type $job_type $extract_raw_manual_dir $script_type jobs1 steps1 j
            multi_jobs $pbs_run_loc $pbs_job_prefix $nci_project $raw_mem $raw_ncpus $queue1 $script $depend_job $depend_type $job_type $extract_raw_manual_dir $script_type jobs2 steps2 jobs1
        } < $list

        # Error collation
        cd $extract_raw_batch_dir
        echo ""
        echo "Preparing error collation for raw data extraction ..."
        depend_job=`sed s/.r-man2// "all_"$pbs_job_prefix"job_ids"`
        depend_type=afterany
        job_type=1
        pbs_job_prefix1=raw_err
        err_type=1
        script=collate_nci_errors.bash
        {
            single_job $pbs_run_loc $pbs_job_prefix1 $nci_project $extract_raw_batch_dir $err_walltime $err_mem $err_ncpus $exp_queue $depend_job $depend_type $job_type $err_type $script
        }

        # clean up PBS job dir
        cd $extract_raw_batch_dir
        rm -rf list* $pbs_job_prefix*"_job_id"
   fi

elif [ $do_raw == no ]; then
    echo "Option to extract raw data not selected."
    echo ""
else
    :
fi


