#!/bin/bash

surfID=191022NHJ
SN=13NHJ
unresampled_data=( CrM1.nii )

SURF_DIR=/sas2/PECON/freesurfer/subjects
LOC_DIR=/sas2/PECON/HJY/CrM/Exp2/${SN}/localizerAnalysis/
MNI_DIR=/sas2/PECON/HJY/CrM/MNItemplate
ROI_DIR=${MNI_DIR}/A1_template
tete=RAI_MNI152_T1_2mm+tlrc #template 

cd ${LOC_DIR}

run=( 1 2 )
name=audLoc
a_seNm=CIN #same session with audLoc

for r in "${run[@]}"
do
@Align_Centers -no_cp -base Inplane_${a_seNm}@CrMvr2+orig. -dset audLoc${r}.nii
done

#preprocessing
afni_proc.py -subj_id ${name} \
-dsets audLoc?.nii \
-volreg_base_dset ${name}1.nii['2'] \
-copy_anat Inplane_${a_seNm}@CrMvr2+orig. \
-regress_censor_motion 0.5 \
-volreg_align_e2a \
-align_opts_aea -cost lpc+ZZ -giant_move \
-regress_censor_outliers 0.1 \
-blocks align volreg mask regress
tcsh -xf proc.${name}

cd ${name}.results
afni & #check alignment 

#-------------- Remove background ------------------# 
run=( 01 02 )
for r in "${run[@]}"
do
3dClipLevel pb01.${name}.r${r}.volreg+orig.HEAD >> clip.txt
3dTstat -mean -prefix r.${r}.base pb01.${name}.r${r}.volreg+orig.HEAD'[0..$]'
done

more clip.txt # check the smallest clip value across all runs
clip=$(cut -f1 -d"," clip.txt | sort -n | head -1) # assign the clip value

#-------------- pb02 Scaling ------------------# 
for r in "${run[@]}"
do
3dcalc -a pb01.${name}.r${r}.volreg+orig. -b r.${r}.base+orig. \
	-expr "(100 * a/b) * step(b-$clip)" -prefix pb02.${name}.r${r}.scaled #remove the space and scale 
done

#-------------- pb03 Detrend & 04 Highpass filter & pb05 Add mean  ------------------# 

for i in "${run[@]}"
do
#detrend & highpass filter
3dDetrend -polort 1 -prefix pb03.${name}.r$i.sc_dt+orig pb02.${name}.r$i.scaled+orig
3dBandpass -prefix pb04.${name}.r${i}.sc_dt_hp 0.01 99999 pb03.${name}.r${i}.sc_dt+orig
#add mean of scaled data (pb02.~)
3dTstat -mean -prefix r.${i}.sc_base pb02.${name}.r${i}.scaled+orig.HEAD'[0..$]'
3dcalc  -a pb04.${name}.r${i}.sc_dt_hp+orig.HEAD -b r.${i}.sc_base+orig.HEAD \
	-expr 'a+b' -prefix pb05.${name}.r${i}.sc_dt_hp_am
done

#-------------- pb06 smoothing  ------------------# 
for i in "${run[@]}"
do
3dmerge -1blur_fwhm 5 -doall -prefix rm.pb06.${name}.r${i}.sc_dt_hp_am_blur    \
pb05.${name}.r${i}.sc_dt_hp_am+orig

# and set boundaries using anat mask 
3dcalc -a rm.pb06.${name}.r${i}.sc_dt_hp_am_blur+orig -b full_mask.${name}+orig. \
-expr 'a*b' -prefix pb06.${name}.r${i}.sc_dt_hp_am_blur+orig

done

rm -f rm.pb06*

#-------------- Regress ------------------# 
3dDeconvolve -input pb06.${name}.r*.sc_dt_hp_am_blur+orig.HEAD            \
    -censor censor_${name}_combined_2.1D                                    \
    -polort 3 							            \
    -num_stimts 8 			                                    \
    -stim_file 1 motion_demean.1D'[0]' -stim_base 1 -stim_label 1 roll      \
    -stim_file 2 motion_demean.1D'[1]' -stim_base 2 -stim_label 2 pitch     \
    -stim_file 3 motion_demean.1D'[2]' -stim_base 3 -stim_label 3 yaw       \
    -stim_file 4 motion_demean.1D'[3]' -stim_base 4 -stim_label 4 dS        \
    -stim_file 5 motion_demean.1D'[4]' -stim_base 5 -stim_label 5 dL     \
    -stim_file 6 motion_demean.1D'[5]' -stim_base 6 -stim_label 6 dP     \
    -stim_times 7 /sas2/PECON/HJY/CrM/localizerOnset/+stim-fixation.txt 'BLOCK(8,1)'           \
    -stim_label 7 Audio                                                    \
    -stim_times 8 /sas2/PECON/HJY/CrM/localizerOnset/+fixation-stim.txt 'BLOCK(8,1)'           \
    -stim_label 8 Silence                                                 \
    -local_times                                                            \
    -gltsym 'SYM: +Audio -Silence'                     \
    -glt_label 1 audio_vs_silence                      \
    -gltsym 'SYM: +Silence -Audio'                     \
    -glt_label 2 silence_vs_audio                      \
    -float                                                                  \
    -jobs 8                                                                 \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                 \
    -x1D_uncensored X.nocensore54.xmat.1D                                      \
    -bucket stats.${name}

#---------------Thresholding stat results ------------------# 
#cp ../AVol@CrMvr2+orig.* .

afni -niml &
suma -spec ${SURF_DIR}/${surfID}/SUMA/${surfID}_both.spec -sv ../SVol@CrMvr2+orig.HEAD

# click on suma window and press 't'
# threshold q=.001

# save contrast as separate file
3dbucket stats.${name}+orig.HEAD[8] -prefix stat.audLocResults

# find thresholding t-value in AFNI and make threshold mask
3dcalc -a stat.audLocResults+orig. -expr "ispositive(a-4.297)" -prefix audLoc_thresmask

####### Align Auditory Masks in MNI space to Native space#######
 
#-----------------------method 1. FLIRT (FNL) -------------------------
#Use Flirt instead...
#ref: thanks to hyssong
#/group_hpc/WMShimLab/2FilmStudy/script/a_preprocessing/preprocessing_edit3.m

mkdir flirt@MNI.results
cd flirt@MNI.results 
flirt_tete=RAI_MNI152_T1_2mm.nii.gz
cp ../../AVol@CrMvr2+orig.* .
#cp ../../rawdata/session1/"${unresampled_data[0]}" .; ../../rawdata/session1/"${unresampled_data[1]}" .
cp ${MNI_DIR}/RAI_orientation/${flirt_tete} .

#AFNI to NIFTI (AVol)
3dAFNItoNIFTI -prefix AVol@CrMvr2.nii AVol@CrMvr2+orig.

flirt -ref ${flirt_tete} -in AVol@CrMvr2.nii.gz -out AVol@MNI.nii -omat normMatrix.mat -dof 12 -cost mutualinfo #normcorr, normmi, leastsq, 
convert_xfm -inverse -omat inv.normMatrix.mat normMatrix.mat

## if you have an error ... 
#!/bin/bash

# Read from specified file, or from standard input
infile="${1:-/sas2/PECON/HJY/CrM/Exp2/13NHJ/localizerAnalysis/audLoc.results/flirt@MNI.results/inv.normMatrix.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile
###

#MNI to Native (rawdata orientation)
ROIs=( TE10_lh TE10_rh TE11_lh TE11_rh TE12_lh TE12_rh ) 
for i in "${ROIs[@]}"
do
flirt -ref AVol@CrMvr2.nii -in ${MNI_DIR}/RAI_orientation/RAI_${i}.nii -out tmp_Native_${i} -applyxfm -init inv.normMatrix.mat
3dcopy tmp_Native_${i}.nii.gz Native_${i} -overwrite 
done

rm -f tmp_Native_*.nii.gz


afni -niml &
suma -spec ${SURF_DIR}/${surfID}/SUMA/${surfID}_both.spec -sv ../../SVol@CrMvr2+orig.HEAD

rm -f ${flirt_tete}
##---------------Method 2. AFNI: Get Inverted Matrix (MNI -> Native)--------------------------# 

##Align to MNI space first 
cd ../
cp ../mainPrep/"${unresampled_data[0]}" .; cp ../mainPrep/"${unresampled_data[1]}" .; 

afni_proc.py -subj_id afni@MNI \
-dsets CrM*.nii \
-blocks align tlrc volreg blur mask regress \
-align_opts_aea -giant_move               \
-tlrc_base ${MNI_DIR}/RAI_orientation/${tete} \
-volreg_base_dset CrM1.nii'[0]' \
-volreg_tlrc_warp \
-copy_anat AVol@CrM+orig. \
-regress_censor_motion 0.5                \
-regress_censor_outliers 0.1

tcsh -xef proc.afni@MNI |& tee output.proc.afni@MNI
cd ./afni@MNI.results/
mv warp.anat.Xat.1D normMatrix.1D
cat_matvec normMatrix.1D -I > inv.normMatrix.1D #make inv. matrix.1D 
cd ../

###or try this....
@auto_tlrc -base MNI_avg152T1+tlrc -input AVol@CrM+orig. -init_xform AUTO_CENTER -no_ss 
###or try this....
align_epi_anat.py -dset1 RAI_MNI152_T1_2mm+tlrc  -dset2 AVol@CrM+orig. -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move

#Convert base dataset format to nii.gz or nii  
3dAFNItoNIFTI -prefix CrM1.nii.gz CrM1.nii 

#MNI to Native (rawdata orientation)
ROIs=( TE10_lh TE10_rh TE11_lh TE11_rh TE12_lh TE12_rh ) 
for i in "${ROIs[@]}"
do
3dAllineate -base normEPI_CrM1.nii.gz -input ${ROI_DIR}/${i}_bin.nii.gz -1Dmatrix_apply ./normMatrix.results/inv.normMatrix.1D \
-mast_dxyz 2 -prefix Native_${i} -overwrite
done

#---------------Remove overlapping voxels from TE Masks  ------------------# 
#combine lh&rh / threshold aud ROIs from MNI once again (I don't know why)

ROIs=( TE10 TE11 TE12 )

for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r Native_${r}_lh+orig. Native_${r}_rh+orig.
3dcalc -a ${r}+orig. -expr 'ispositive(a-0.5)' -prefix ${r}_bin
done

# exclude overlapping voxels
# count voxels of each ROI: TE1.2 has the smallest number of voxels
# 3dBrickStat -sum {ROIname}

# TE1.2 & TE1.1 -> exclude from TE1.1 mask
mv TE11_bin+orig.HEAD tmp.TE11+orig.HEAD; mv TE11_bin+orig.BRIK tmp.TE11+orig.BRIK
3dcalc -a tmp.TE11+orig -b TE12_bin+orig -expr 'a+b' -prefix tmp.TE11+12+orig
3dcalc -a tmp.TE11+12+orig -expr 'ispositive(a-1.5)' -prefix mask.TE11+12+orig
3dcalc -a tmp.TE11+orig -b mask.TE11+12+orig -expr 'ispositive(a-b)' \
	-prefix TE11_bin+orig #from 1029 to 990

# TE1.1 & TE1.0 -> exclude from TE1.0 mask
mv TE10_bin+orig.HEAD tmp.TE10+orig.HEAD; mv TE10_bin+orig.BRIK tmp.TE10+orig.BRIK
3dcalc -a tmp.TE10+orig -b TE11_bin+orig -expr 'a+b' -prefix tmp.TE11+10+orig
3dcalc -a tmp.TE11+10+orig -expr 'ispositive(a-1.5)' -prefix mask.TE11+10+orig
3dcalc -a tmp.TE10+orig -b mask.TE11+10+orig -expr 'ispositive(a-b)' \
	-prefix TE10_bin+orig #from 1195 to 766

#Count voxels again 
for r in "${ROIs[@]}"
do
3dBrickStat -sum ${r}_bin+orig.
3dresample -dxyz 2 2 2 -inset ${r}_bin+orig. -prefix ${r}_bin+orig. -overwrite 
done 

rm -f tmp* mask.TE*

cd ../
#---------------ThresMask X Native_A1_Template (Final Masks) ------------------# 

# threshold with functional contrast map
for r in "${ROIs[@]}"
do
3dcalc -a flirt@MNI.results/${r}_bin+orig.HEAD -b ./audLoc_thresmask+orig -expr 'ispositive(a*b)' -prefix ${r}_masked_tmp 
3dcalc -a ${r}_masked_tmp+orig.HEAD -b ./full_mask.${name}+orig.HEAD \
-expr 'ispositive(a*b)' -prefix  ${r}_fmasked.nii
done

# Merge all Masks -> A1! 
3dmerge -gmax -prefix A1_fmasked.nii TE10_fmasked.nii TE11_fmasked.nii TE12_fmasked.nii 

rm -f *_masked_tmp*
mv *_fmasked.nii ../../ROImasks/




