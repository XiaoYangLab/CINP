#!/bin/bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw


for subj in `ls $subjDIR`; do

mkdir -p $subjDIR/${subj}/mytbss

###################cp meas
for i in FA MD MO L1 L2 L3 ; do
cp $subjDIR/${subj}/DTI_preproc/dti_fit_${i}.nii.gz $subjDIR/${subj}/mytbss/dti_fit_${i}.nii.gz
done

###################tbss for FA
cd $subjDIR/${subj}/mytbss
tbss_1_preproc dti_fit_FA.nii.gz
tbss_2_reg -T
tbss_3_postreg -T
tbss_4_prestats 0.2
${FSLDIR}/bin/fslstats -K ${FSLDIR}/data/atlases/JHU/JHU-ICBM-labels-1mm $subjDIR/${subj}/mytbss/stats/all_FA_skeletonised.nii.gz -M >$subjDIR/${subj}/mytbss/stats/JHUrois_FA.txt

###################tbss for MD MO L1 L2 L3
for i in MD MO L1 L2 L3 ; do
${FSLDIR}/bin/applywarp --rel -i $subjDIR/${subj}/mytbss/dti_fit_${i} -o $subjDIR/${subj}/mytbss/stats/all_${i} -r ${FSLDIR}/data/standard/FMRIB58_FA_1mm -w $subjDIR/${subj}/mytbss/FA/dti_fit_FA_FA_to_target_warp
${FSLDIR}/bin/fslmaths $subjDIR/${subj}/mytbss/stats/all_${i} -mas $subjDIR/${subj}/mytbss/stats/mean_FA_skeleton_mask $subjDIR/${subj}/mytbss/stats/all_${i}_skeletonised
${FSLDIR}/bin/fslstats -K ${FSLDIR}/data/atlases/JHU/JHU-ICBM-labels-1mm $subjDIR/${subj}/mytbss/stats/all_${i}_skeletonised.nii.gz -M >$subjDIR/${subj}/mytbss/stats/JHUrois_${i}.txt
done

done
