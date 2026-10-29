#!/usr/bin/env bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw
#subjDIR='/media/kame/Elements SE/Single_APDs/Risperidone'
#subjDIR='/media/kame/Elements SE/Single_APDs/Olanzapine'
subjDIR='/media/kame/Elements SE/Single_APDs/Aripiprazole'

for subj in `ls "$subjDIR"`
do
mkdir -p "$subjDIR"/${subj}/DTI_nii
dcm2niix -o "$subjDIR"/${subj}/DTI_nii "$subjDIR"/${subj}/DTI

done
