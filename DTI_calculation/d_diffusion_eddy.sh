#!/bin/bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
#error_subj=("SZ_01_0021" "SZ_01_0022")
#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw
#error_subj=("SZ_10_0028" "SZ_10_0039")

subjDIR=/media/kame/pku6/Single_APDs/Risperidone
error_subj=("1251691")
#subjDIR=/media/kame/pku6/Single_APDs/Olanzapine
#subjDIR=/media/kame/pku6/Single_APDs/Aripiprazole

#template=""
#tool=""
#parameter=""

for subj in `ls "$subjDIR"`
#for subj in "${error_subj[@]}"
do
mkdir -p "$subjDIR"/${subj}/DTI_preproc
cd "$subjDIR"/${subj}/DTI_preproc

# prepare data for eddy, dtifit and bedpostx.
cp "$subjDIR"/${subj}/DTI_nii/*.bval ./AP.bval
cp "$subjDIR"/${subj}/DTI_nii/*.bvec ./AP.bvec
cp "$subjDIR"/${subj}/DTI_nii/*.nii ./AP.nii
cp "$subjDIR"/${subj}/DTI_nii/*.json ./AP.json


# extract b0 image
${FSLDIR}/bin/fslroi AP.nii b0 0 1

# extract bo mask
${FSLDIR}/bin/bet2 b0.nii.gz b0_brain -m -f 0.2 

########
printf "0 1 0 0.05\n0 -1 0 0.05" >> ./acqparams.txt
indx=""
n=`${FSLDIR}/bin/fslval AP.nii dim4`
for ((i=1;i<=${n};i++));do
indx="$indx 1"
done
echo $indx > ./eddy_index.txt

# eddy
${FSLDIR}/bin/eddy_cuda10.2 --imain=./AP.nii --mask=./b0_brain_mask.nii.gz --acqp=./acqparams.txt --index=./eddy_index.txt --out=./AP_eddy --bvecs=./AP.bvec --bvals=./AP.bval

# rotate by eddy correction
#fdt_rotate_bvecs AP.bvec AP_rotate.bvec AP_eddy.ecclog

# ditfit for FA, MD, MO, L1, L2, L3
${FSLDIR}/bin/dtifit -k ./AP_eddy.nii.gz -m ./b0_brain_mask.nii.gz -r ./AP_eddy.eddy_rotated_bvecs -b ./AP.bval -o dti_fit --sse 


#####################bb_get_b0 from ukbiobank pipeline
## -n is the number of b0 in dwi
# python $tool/bb_get_b0s.py -i ./AP.nii -o ./total_B0_AP.nii -n 1 -l 100 -a bvals

#####################choose_bestB0 from ukbiobank pipeline
# $tool/bb_choose_bestB0 ./total_B0_AP.nii ./B0_AP.nii
# bestInd='cat ./total_B0_AP_best_index.txt'

#####################run topup for dMRI
# topupValue='echo 0.67*127*0.001|bc' #echo spacing * epi facter length -1 * 0.001
# printf '0 -1 0 $topupValue' >> ./acqparams.txt
# topup --imain=./total_B0_AP.nii --datain=./acqparams.txt --config=$parameter /
# topup_b02b0.cnf --out=./topup_out --fout=./topup_fout --jacout=./topup_jacout -v

## Get corrected B=0
# ${FSLDIR}/bin/applytopup --imain=./B0_AP.nii --datain=./acqparams.txt --inindex=1 --
# topup=./topup_out --out=./B0_AP_corr_tmp --method=jac
# python $tool/bb_mask_negatives_4D.py -i ./B0_AP_corr_tmp.nii -o ./B0_AP_mask.nii -z
# True

# ${FSLDIR}/bin/fslmaths ./B0_AP_corr_tmp.nii -mul ./B0_AP_mask.nii ./AP_masked_tmp.nii
# ${FSLDIR}/bin/fslmaths ./B0_AP_corr_tmp.nii -add ./AP_masked_tmp.nii ./B0_AP_fixed_tmp.nii
# mv ./B0_AP_fixed_tmp.nii ./topup_iout.nii

# $tool/bb_GDC --workingdir = ./topup_iout_GDC/ --in=./topup_iout.nii --out=./topup_iout_ud.nii.gz --owarp=./topup_iout_ud_warp.nii
# $tool/bb_epi_reg --epi=./B0_AP.nii --t1=../T1_nii/*.nii --t1brain=../T1_nii/T1_brain.nii --out=./B0_to_T1 
# --wmseg=../T1/T1fast/T1_brain_WM_mask.nii

## generate a mask for topup output by inverting the previous registration and applying it to T1 brain mask
# ${FSLDIR}/bin/convert_xfm -omat ./T1_to_B0.mat -inverse ./B0_to_T1.mat
# ${FSLDIR}/bin/fslmaths ../T1/T1fast/T1_brain_mask.nii -thr 0.1 -kernel sphere 1.1 -dilF -bin 
# -fillh ../T1/T1fast/T1_brain_mask_dil.nii
# ${FSLDIR}/bin/flirt -in ../T1/T1fast/T1_brain_mask_dil.nii -ref ./topup_iout.nii -applyxfm -init ./T1_to_B0.mat 
# -out ./topup_mask_ud.nii -interp trilinear
# ${FSLDIR}/bin/fslmaths ./topup_mask_ud.nii -thr 0.25 -bin -fillh ./topup_mask_ud.nii

## warp the dilated T1 brain mask to the gradient distored space by inverting the mat
# ${FSLDIR}/bin/invwarp --ref=./topup_iout_ud.nii -w ./topup_iout_ud_warp.nii -o ./topup_iout_mean_ud_inv_warp.nii
# ${FSLDIR}/bin/applywarp --rel -i ../T1/T1fast/T1_brain_mask_dil.nii -r ./topup_iout.nii -o ./topup_mask.nii 
# --premat=./T1_to_B0.mat -w ./topup_iout_mean_ud_inv_warp.nii --interp=trilinear
# ${FSLDIR}/bin/fslmaths ./topup_mask.nii -thr 0.25 -bin -fillh ./topup_mask.nii

## eddy
# ${FSLDIR}/bin/eddy --imain=./AP.nii --mask=./topup_mask.nii --topup=./topup_out --acqp=./acqparams.txt --index=./eddy_index.txt --#bvecs=./bvecs --bvals=./bvals --out=./data --ref_scan_no=$bestInd --flm=quadratic --resamp=jac --slm=linear --fwhm=2 --ff=5 --#sep_offs_move --nvoxhp=1000 --very_verbose --repol --rms

#dtifit for FA MD MO
#dtifit -k ./data_ud_1_shell -m ./topup_mask_ud -r ./data_ud_1_shell.bvec -b ./data_ud_1_shell.bval -o ./dti


done
