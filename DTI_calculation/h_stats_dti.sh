#!/bin/bash

#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_PKU6_MRI_baseline_raw
#resultsDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/results/PKU6
#subjDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/0_dataset/APDs/SZ_ZMD_MRI_baseline_raw
#resultsDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/results/ZMD

#subjDIR=/media/kame/pku6/Single_APDs/Risperidone
#resultsDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/results/Risperidone
subjDIR=/media/kame/pku6/Single_APDs/Olanzapine
resultsDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/results/Olanzapine
#subjDIR=/media/kame/pku6/Single_APDs/Aripiprazole
#resultsDIR=/media/kame/Data/_workings/projects/SZImagingGenetic/results/Aripiprazole

#####################setting the string of NaN in case there is a problem
numVars="27"
nanResults=$(printf "NaN %.0s" $(seq 1 $numVars))

#####################extract for each subject
for subj in `ls $subjDIR`
do
echo "Processing subject: ${subj}"
#####################check all tracts
ntracts=`echo $subjDIR/${subj}/autoptx_preproc/tracts/*/tracts/tractsNorm.nii.gz | wc -w`
if [ $ntracts != 27 ]; then
>&2 echo only found $ntracts tracts
fi

#####################hack to ameliorate partial misaligment of mcp
${FSLDIR}/bin/flirt -in ${FSLDIR}/data/standard/FMRIB58_FA_1mm -ref $subjDIR/${subj}/DTI_preproc/refVol -applyisoxfm 1 -out $subjDIR/${subj}/autoptx_preproc/FMRIB58_FA_1mm_temp
fslmaths $subjDIR/${subj}/autoptx_preproc/FMRIB58_FA_1mm_temp -thr 1000 -bin -ero -ero -mul $subjDIR/${subj}/autoptx_preproc/tracts/mcp/tracts/tractsNorm $subjDIR/${subj}/autoptx_preproc/tracts/mcp/tracts/tractsNorm -odt float
/bin/rm -f $subjDIR/${subj}/autoptx_preproc/FMRIB58_FA_1mm_temp.nii.gz

#####################normalize by 99% threshold
for p in `ls $subjDIR/${subj}/autoptx_preproc/tracts`; do
if [ -f $subjDIR/${subj}/autoptx_preproc/tracts/$p/tracts/tractsNorm.nii.gz ]; then
t=`fslstats $subjDIR/${subj}/autoptx_preproc/tracts/$p/tracts/tractsNorm -P 99`
fslmaths $subjDIR/${subj}/autoptx_preproc/tracts/$p/tracts/tractsNorm -min $t -mul 1e6 -div $t -thr 1e5 $subjDIR/${subj}/autoptx_preproc/tmp_$p
fi
done

#####################merge all tractNorm
pp=$(ls $subjDIR/${subj}/autoptx_preproc/tmp_*.nii.gz)
fslmerge -t $subjDIR/${subj}/autoptx_preproc/tmp_all $pp
tractnorm=`fslstats -t $subjDIR/${subj}/autoptx_preproc/tmp_all -m`

#####################launch statistics for autoptx
tractnames=($(ls -d $subjDIR/${subj}/autoptx_preproc/tracts/*/ | xargs -n 1 basename))
for i in FA MD MO L1 L2 L3 ; do
/bin/rm -f $subjDIR/${subj}/autoptx_preproc/APTX_${i}.txt

${FSLDIR}/bin/flirt -in $subjDIR/${subj}/DTI_preproc/dti_fit_${i} -ref $subjDIR/${subj}/DTI_preproc/dti_fit_${i} -applyisoxfm 1 -out $subjDIR/${subj}/autoptx_preproc/tmp_all_${i}  #$subjDIR/${subj}/DTI_preproc/refVol
fslmaths $subjDIR/${subj}/autoptx_preproc/tmp_all -mul $subjDIR/${subj}/autoptx_preproc/tmp_all_${i} $subjDIR/${subj}/autoptx_preproc/tmp
tractmeas=`fslstats -t $subjDIR/${subj}/autoptx_preproc/tmp -m | sed 's/-/_/g'`
jmax=`fslnvols $subjDIR/${subj}/autoptx_preproc/tmp_all`
for (( j=1; j <= $jmax; ++j )) ; do
tm=`echo $tractmeas | cut -d " " -f$j`
tn=`echo $tractnorm | cut -d " " -f$j`
fibername=${tractnames[$((j-1))]}
echo "${fibername} $(echo "10 k $tm $tn / p" | dc -)" >> $subjDIR/${subj}/autoptx_preproc/APTX_${i}.txt
done

done
/bin/rm $subjDIR/${subj}/autoptx_preproc/tmp*.nii.gz

done

########################
mkdir -p $resultsDIR
mkdir -p $resultsDIR/fibers

for subj in `ls $subjDIR`
do
echo -e "ID\tFiber\tFA\tMD\tMO\tL1\tL2\tL3" > "$resultsDIR/fibers/${subj}_fibers.csv"

for fiber in $(ls -d $subjDIR/${subj}/autoptx_preproc/tracts/*/ | xargs -n 1 basename); do
fa=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_FA.txt")
md=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_MD.txt")
mo=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_MO.txt")
l1=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_L1.txt")
l2=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_L2.txt")
l3=$(awk -v fiber="$fiber" '$1 == fiber {print $2}' "$subjDIR/${subj}/autoptx_preproc/APTX_L3.txt")

echo -e "${subj}\t${fiber}\t$fa\t$md\t$mo\t$l1\t$l2\t$l3" >> "$resultsDIR/fibers/${subj}_fibers.csv"
done
done








