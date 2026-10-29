#!/bin/bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
#error_subj=("SZ_01_0021" "SZ_01_0022")
#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw
#error_subj=("SZ_10_0028" "SZ_10_0039")

subjDIR=/media/kame/pku6/Single_APDs/Risperidone
error_subj=("1251691")
#subjDIR=/media/kame/pku6/Single_APDs/Olanzapine
#subjDIR=/media/kame/pku6/Single_APDs/Aripiprazole

for subj in `ls "$subjDIR"`
#for subj in "${error_subj[@]}"
do

${FSLDIR}/bin/immv "$subjDIR"/${subj}/DTI_preproc/AP_eddy.nii.gz "$subjDIR"/${subj}/DTI_preproc/data
${FSLDIR}/bin/immv "$subjDIR"/${subj}/DTI_preproc/b0_brain_mask "$subjDIR"/${subj}/DTI_preproc/nodif_brain_mask
\mv "$subjDIR"/${subj}/DTI_preproc/AP.bvec "$subjDIR"/${subj}/DTI_preproc/bvecs
#\mv "$subjDIR"/${subj}/DTI_preproc/AP_eddy.eddy_rotated_bvecs "$subjDIR"/${subj}/DTI_preproc/bvecs
\mv "$subjDIR"/${subj}/DTI_preproc/AP.bval "$subjDIR"/${subj}/DTI_preproc/bvals

##################### estimate registration to standard space
${FSLDIR}/bin/fsl_reg "$subjDIR"/${subj}/DTI_preproc/dti_fit_FA ${FSLDIR}/data/standard/FMRIB58_FA_1mm "$subjDIR"/${subj}/DTI_preproc/nat2std -e -FA
${FSLDIR}/bin/invwarp -r "$subjDIR"/${subj}/DTI_preproc/dti_fit_FA -w "$subjDIR"/${subj}/DTI_preproc/nat2std_warp -o "$subjDIR"/${subj}/DTI_preproc/std2nat_warp

##################### launch bedpostx
${FSLDIR}/bin/bedpostx_gpu "$subjDIR"/${subj}/DTI_preproc

##################### create native space reference volume in 1mm cubic resolution
${FSLDIR}/bin/flirt -in "$subjDIR"/${subj}/DTI_preproc/dti_fit_FA -ref "$subjDIR"/${subj}/DTI_preproc/dti_fit_FA -applyisoxfm 1 -out "$subjDIR"/${subj}/DTI_preproc/refVol

${FSLDIR}/bin/fslmaths "$subjDIR"/${subj}/DTI_preproc/refVol -mul 0 "$subjDIR"/${subj}/DTI_preproc/refVol



#############################Making bedpostx directory structure

#mkdir -p $subjDIR/${subj}/DTI_preproc/bedpostx
#mkdir -p $subjDIR/${subj}/DTI_preproc/bedpostx/diff_parts
#mkdir -p $subjDIR/${subj}/DTI_preproc/bedpostx/logs
#mkdir -p $subjDIR/${subj}/DTI_preproc/bedpostx/logs/logs_gpu
#mkdir -p $subjDIR/${subj}/DTI_preproc/bedpostx/xfms

#cp $subjDIR/${subj}/DTI_preproc/AP.bvec $subjDIR/${subj}/DTI_preproc/AP.bval $subjDIR/${subj}/DTI_preproc/bedpostx 

#${FSLDIR}/bin/imcp $subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz $subjDIR/${subj}/DTI_preproc/bedpostx

#############################xfibres gpu version
#${FSLDIR}/bin/split_parts_gpu $subjDIR/${subj}/DTI_preproc/AP_eddy.nii $subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz $subjDIR/${subj}/DTI_preproc/AP.bval $subjDIR/${subj}/DTI_preproc/AP.bvec NULL 0 1 $subjDIR/${subj}/DTI_preproc/bedpostx

#nvox=`${FSLDIR}/bin/fslstats $subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz -V | cut -d ' ' -f1`
#echo $nvox > $subjDIR/${subj}/DTI_preproc/bedpostx/nvox.txt

#############################bedpostx
#opts="--nf=3 --fudge=1 --bi=1000 --nj=1250 --se=25 --model=2 --cnonlinear"
#nvox=`cat $subjDIR/${subj}/DTI_preproc/bedpostx/nvox.txt`

#${FSLDIR}/bin/xfibres_gpu --data=$subjDIR/${subj}/DTI_preproc/AP_eddy.nii --mask=$subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz -b $subjDIR/${subj}/DTI_preproc/AP.bval -r $subjDIR/${subj}/DTI_preproc/AP.bvec --forcedir --logdir=$subjDIR/${subj}/DTI_preproc/bedpostx/diff_parts/data_part_0000 $opts $subjDIR/${subj}/DTI_preproc 0 1 $nvox

#opts="--nf=3 --fudge=1 --bi=1000 --nj=1250 --se=25 --model=2 --cnonlinear"
#nvox=`cat $subjDIR/${subj}/DTI_preproc/bedpostx/nvox.txt`
#mv $subjDIR/${subj}/DTI_preproc/bedpostx $subjDIR/${subj}/DTI_preproc.bedpostX
#cp $subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz $subjDIR/${subj}/DTI_preproc.bedpostX/nodif_brain_mask.nii.gz
#cp $subjDIR/${subj}/DTI_preproc/b0_brain_mask.nii.gz $subjDIR/${subj}/DTI_preproc/nodif_brain_mask.nii.gz

#${FSLDIR}/bin/bedpostx_postproc_gpu.sh --data=$subjDIR/${subj}/DTI_preproc/AP_eddy.nii --mask=$subjDIR/${subj}/DTI_preproc.bedpostX/nodif_brain_mask.nii.gz -b $subjDIR/${subj}/DTI_preproc.bedpostX/AP.bval -r $subjDIR/${subj}/DTI_preproc.bedpostX/AP.bvec --forcedir --logdir=$subjDIR/${subj}/DTI_preproc.bedpostX/diff_parts $opts $nvox 1 $subjDIR/${subj}/DTI_preproc "${FSLDIR}"

done

echo "all BedpostX jobs have finished."
