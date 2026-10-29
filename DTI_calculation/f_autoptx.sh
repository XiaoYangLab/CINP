#!/bin/bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
#error_subj=("" "")
#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw
#error_subj=("SZ_10_0028" "SZ_10_0039")

#subjDIR=/media/kame/pku6/Single_APDs/Risperidone
#subjDIR=/media/kame/pku6/Single_APDs/Olanzapine
subjDIR=/media/kame/pku6/Single_APDs/Aripiprazole

structures=/media/kame/Data/_workings/projects/SZImagingGenetic/3_preprocessing/_autoptx/structureList

#####################launch Tractography

while read structstring; do
    struct=`echo $structstring | awk '{print $1}'`
    nseed=`echo $structstring | awk '{print $2}'`
    wallt=`echo $structstring | awk '{print $3}' | tr ':' ' '`
    hours=`echo $wallt | awk '{print $1}'`
    minutes=`echo $wallt | awk '{print $2}'`
    time=$(( 60 * hours + minutes ))
    masks=/media/kame/Data/_workings/projects/SZImagingGenetic/3_preprocessing/_autoptx/protocols/$struct

    ################ is there a stop criterion defined in the protocol for this struct?
    if [ -e $masks/stop.nii.gz ]; then
    useStop=1
    else
    useStop=0
    fi

    ################ apply seed multiplier if set
    nSeed=1000
    if [ "$nseed" != "" ]; then
    nSeed=$(echo "scale=0; $nSeed * $nseed / 1"|bc)
    fi

    ############### run subject list
    for subj in `ls $subjDIR`
    #for subj in "${error_subj[@]}"
    do

    echo "running automated tractography of subject ${subj}, structure ${struct}, using $nSeed seeds per voxel."
    mkdir -p $subjDIR/${subj}/autoptx_preproc
    mkdir -p $subjDIR/${subj}/autoptx_preproc/logs
    mkdir -p $subjDIR/${subj}/autoptx_preproc/tracts
    mkdir -p $subjDIR/${subj}/autoptx_preproc/tracts/$struct
    
    ################# cleanup possible previous run
    \rm -f $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/waytotal
    imrm $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/tractsNorm
    \rm -rf $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts
    
    ################ does the protocol defines a second run with inverted seed / target masks
    if [ -e $masks/invert ]; then
    symtrack=1
    rm -f $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv/waytotal
    else
    symtrack=0
    fi

    ################# transform masks
    ${FSLDIR}/bin/applywarp -r $subjDIR/${subj}/DTI_preproc/refVol -w $subjDIR/${subj}/DTI_preproc/std2nat_warp -i $masks/seed -o $subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed -d float
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed -thr 0.1 -bin $subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed -odt char

    ${FSLDIR}/bin/applywarp -r $subjDIR/${subj}/DTI_preproc/refVol -w $subjDIR/${subj}/DTI_preproc/std2nat_warp -i $masks/target -o $subjDIR/${subj}/autoptx_preproc/tracts/$struct/target -d float
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/target -thr 0.1 -bin $subjDIR/${subj}/autoptx_preproc/tracts/$struct/target -odt char

    ${FSLDIR}/bin/applywarp -r $subjDIR/${subj}/DTI_preproc/refVol -w $subjDIR/${subj}/DTI_preproc/std2nat_warp -i $masks/exclude -o $subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -d float
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -thr 0.1 -bin $subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -odt char

    if [ "$useStop" -eq "1" ]; then
    ${FSLDIR}/bin/applywarp -r $subjDIR/${subj}/DTI_preproc/refVol -w $subjDIR/${subj}/DTI_preproc/std2nat_warp -i $masks/stop -o $subjDIR/${subj}/autoptx_preproc/tracts/$struct/stop -d float
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/stop -thr 0.1 -bin $subjDIR/${subj}/autoptx_preproc/tracts/$struct/stop -odt char
    fi


    ################### process structure
    if [ "$useStop" -eq "1" ]; then
    ${FSLDIR}/bin/probtrackx2_gpu -s $subjDIR/${subj}/DTI_preproc.bedpostX/merged -m $subjDIR/${subj}/DTI_preproc/nodif_brain_mask -x $subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed -o density --waypoints=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/target --stop=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/stop --nsamples=${nSeed} --opd --dir=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts --avoid=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -l --forcedir --cthr=0.1 --steplength=1
    if [ "$symtrack" -eq "1" ]; then
    ${FSLDIR}/bin/probtrackx2_gpu -s $subjDIR/${subj}/DTI_preproc.bedpostX/merged -m $subjDIR/${subj}/DTI_preproc/nodif_brain_mask -x $subjDIR/${subj}/autoptx_preproc/tracts/$struct/target -o density --waypoints=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed --stop=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/stop --nsamples=${nSeed} --opd --dir=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv --avoid=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -l --forcedir --cthr=0.1 --steplength=1
    fi

    else
    ${FSLDIR}/bin/probtrackx2_gpu -s $subjDIR/${subj}/DTI_preproc.bedpostX/merged -m $subjDIR/${subj}/DTI_preproc/nodif_brain_mask -x $subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed -o density --waypoints=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/target --nsamples=${nSeed} --opd --dir=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts --avoid=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -l --forcedir --cthr=0.1 --steplength=1
    if [ "$symtrack" -eq "1" ]; then
    ${FSLDIR}/bin/probtrackx2_gpu -s $subjDIR/${subj}/DTI_preproc.bedpostX/merged -m $subjDIR/${subj}/DTI_preproc/nodif_brain_mask -x $subjDIR/${subj}/autoptx_preproc/tracts/$struct/target -o density --waypoints=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/seed --nsamples=${nSeed} --opd --dir=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv --avoid=$subjDIR/${subj}/autoptx_preproc/tracts/$struct/exclude -l --forcedir --cthr=0.1 --steplength=1
    fi

    fi

    ######################### merge runs for forward and inverted tractography runs
    if [ "$symtrack" = "1" ]; then
    ${FSLDIR}/bin/immv $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/density $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv/fwDensity
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv/fwDensity -add $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv/density $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/density
    way1=`cat $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/waytotal | sed 's/e/\\*10^/' | tr -d '+' `
    \rm -f $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/waytotal
    way2=`cat $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tractsInv/waytotal | sed 's/e/\\*10^/' | tr -d '+' `
    way=$(echo "scale=5; $way1 + $way2 "|bc)
    echo $way > $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/waytotal 
    fi

    ########################## perform normalisation for waytotal
    waytotal=`cat $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/waytotal`
    ${FSLDIR}/bin/fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/density -div $waytotal -range $subjDIR/${subj}/autoptx_preproc/tracts/$struct/tracts/tractsNorm -odt float

    echo "containing ${waytotal} valid streamlines\n "	
    done
done < $structures




