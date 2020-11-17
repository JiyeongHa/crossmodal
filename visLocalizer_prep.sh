#!/bin/bash

surfID=190917JJY
SN=14JJY
retID=181
#cpROIs=( lh.V1 lh.V2 lh.V3 rh.V1 rh.V2 rh.V3 lh.fovea2.5dg lh.para+fovea2.5dg rh.fovea2.5dg rh.para+fovea2.5dg )
cpROIs=( lh.V1 lh.V2d lh.V3d lh.V2v lh.V3v rh.V1 rh.V2d rh.V3d rh.V2v rh.V3v lh.fovea3.5dg lh.para+fovea3.5dg rh.fovea3.5dg rh.para+fovea3.5dg lh.FFA lh.PPA lh.MT rh.FFA rh.PPA rh.MT )


SURF_DIR=/sas2/PECON/freesurfer/subjects
LOC_DIR=/sas2/PECON/HJY/CrM/Exp2/${SN}/localizerAnalysis/
RET_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/Img_data
ROI_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/roi/jyh

cd ${LOC_DIR}

run=( 1 2 )
name=visLoc
v_seNm=CIN #same session with visLoc

#visual Localizer
for r in "${run[@]}"
do
@Align_Centers -no_cp -base Inplane_${v_seNm}@CrMvr2+orig. -dset visLoc${r}.nii
done

#preprocessing
afni_proc.py -subj_id ${name} \
-dsets visLoc?.nii \
-copy_anat Inplane_${v_seNm}@CrMvr2+orig. \
-regress_censor_motion 0.5 \
-volreg_align_e2a \
-align_opts_aea -cost lpc+ZZ \
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
    -stim_label 7 Visual                                                   \
    -stim_times 8 /sas2/PECON/HJY/CrM/localizerOnset/+fixation-stim.txt 'BLOCK(8,1)'           \
    -stim_label 8 Fixation                                                \
    -local_times                                                            \
    -gltsym 'SYM: +Visual -Fixation'                     \
    -glt_label 1 visual_vs_fixation                      \
    -gltsym 'SYM: +Fixation -Visual'                     \
    -glt_label 2 fixation_vs_visual                      \
    -float                                                                  \
    -jobs 8                                                                 \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                 \
    -x1D_uncensored X.nocensore54.xmat.1D                                      \
    -bucket stats.${name}

#---------------Thresholding stat results ------------------# 
cp ../AVol@CrMvr2+orig.* .

afni -niml &
suma -spec ${SURF_DIR}/${surfID}/SUMA/${surfID}_both.spec -sv ../SVol@CrMvr2+orig.HEAD

# click on suma window and press 't'
# draw ROI and save as 1D.roi
# threshold q=.001

# save contrast as separate file
3dbucket stats.${name}+orig.HEAD[8] -prefix stat.visLocResults

# find thresholding t-value in AFNI and make threshold mask
3dcalc -a stat.visLocResults+orig. -expr "ispositive(a-4.156)" -prefix visLoc_thresmask

#---------------Copy ROI from Retinotopy data ------------------# 
mkdir roi
cd roi

#cp roi from Retinotopy data or draw on your own
for i in "${cpROIs[@]}"
do
cp ${ROI_DIR}/${i}.1D.roi .
done

#run ROI2Vol.sh lh and rh -> this process aligns ROIs to our basegrid! 
cp /home/jiyeongha/Script/CrMexp2/ROI2Vol.sh .; 
./ROI2Vol.sh lh ${surfID} ${retID}
./ROI2Vol.sh rh ${surfID} ${retID}

ROIs=( V1 )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r lh.${r}+orig. rh.${r}+orig.
done

ROIs=( V2 V3 )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix ${r} lh.${r}d+orig. rh.${r}d+orig. lh.${r}v+orig. rh.${r}v+orig.
done

ROIs=( FFA PPA MT )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r.nii lh.${r}+orig. rh.${r}+orig.
done

# V1 & V2 -> exclude from V1 mask
mv V1+orig.HEAD tmp.V1+orig.HEAD; mv V1+orig.BRIK tmp.V1+orig.BRIK
3dcalc -a tmp.V1+orig -b V2+orig -expr 'a+b' -prefix tmp.V1+V2+orig
3dcalc -a tmp.V1+V2+orig -expr 'ispositive(a-1.5)' -prefix mask.V1+V2+orig
3dcalc -a tmp.V1+orig -b mask.V1+V2+orig -expr 'ispositive(a-b)' \
	-prefix V1+orig

# V2 & V3 -> exclude from V2 mask
mv V2+orig.HEAD tmp.V2+orig.HEAD; mv V2+orig.BRIK tmp.V2+orig.BRIK
3dcalc -a tmp.V2+orig -b V3+orig -expr 'a+b' -prefix tmp.V2+V3+orig
3dcalc -a tmp.V2+V3+orig -expr 'ispositive(a-1.5)' -prefix mask.V2+V3+orig
3dcalc -a tmp.V2+orig -b mask.V2+V3+orig -expr 'ispositive(a-b)' \
	-prefix V2+orig

# V3 & V4(V4v+VO) -> exclude from V3 mask
# 3dmerge -gmax -prefix V4 V4v+orig. VO+orig.
# mv V3+orig.HEAD tmp.V3+orig.HEAD
# mv V3+orig.BRIK tmp.V3+orig.BRIK
# 3dcalc -a tmp.V3+orig -b V4+orig -expr 'a+b' -prefix tmp.V3+V4+orig
# 3dcalc -a tmp.V3+V4+orig -expr 'ispositive(a-1.5)' -prefix mask.V3+V4+orig
# 3dcalc -a tmp.V3+orig -b mask.V3+V4+orig -expr 'ispositive(a-b)' \
#	-prefix V3+orig


rm -f tmp.* mask.V*


#---------------ThresMask X ROI (Final Masks) ------------------# 

ROIs=( V1 V2 V3 )
# threshold with functional contrast map
for r in "${ROIs[@]}"
do
#3dresample -dxyz 2.0 2.0 2.0 -prefix rsmp_${r} -input ${r}+orig. #->match ROI voxel size - thesmask voxel size 
3dcalc -a ${r}+orig.HEAD -b ../${name}_thresmask+orig -expr 'ispositive(a*b)' \
-prefix ${r}_masked #
3dcalc -a ${r}_masked+orig.HEAD -b ../full_mask.${name}+orig.HEAD \
-expr 'ispositive(a*b)' -prefix  ${r}_fmasked.nii
done

rm -f *_masked+orig.*

#when changing roi files to volumetric ones, final result voxels come out of the original space a bit.
#To prevent the mismatch above, multiply volumetric ROI files to final full mask (full mask X ROI )
#This will keep roi files stay in the original boundaries 


##--------------- Make foveal & parafoveal ROIs ------------------##

#change dset to ecc, draw foveal & parafoveal rois here
#change Cmp to afnin16
#My stim was 6.75 degree, and each retino-ecc ring was .61 degree. so, each color represents .61
#so total number of rings I should consider should be at least 11 (11.06)
#In theory, fovea < .6-.75 < parafovea < 5
#36.89 - 1 degree (visual angle)
#final: fovea < 2.5 degree ( 360-73.78=286.22 )
#1 ring: 6.25% , 1 degree = 10.25%
#boundary: 360-249=111


#091719 jyh
#360/16 = 22.5 = 0.61 visual angle
#periphery: 6.75 ->so, at least 11 rings, which correspond to 6.71 -> 22.5*11 = 247.5
#fovea: 3 degree -> at least 5 rings, which correspond to 3.05 -> 22.5*5 = 112.5
#fovea: 3.5 -> 6 rings, 3.66 -> 22.5*6 = 135
#360-247.5 = 112.5
#360 - 135 = 225
# 360 - 22.5 = 337.5
# Experiment1 fixation mask size: 2 
# 360 - 67.5 (1.83 dg) = 292.5



ROIs=( fovea2.5dg para+fovea2.5dg )
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r lh.${r}+orig. rh.${r}+orig.
done

#seperate parafovea mask -> exclude the foveal mask part from parafoveal mask
3dcalc -a para+fovea2.5dg+orig. -b fovea2.5dg+orig. -expr 'a-b' -prefix periphery2.5dg
done

rm -f tmp.* mask.V*

ROIs=( V1 V2 V3 )
# threshold with functional contrast map
for r in "${ROIs[@]}"
do
3dcalc -a ${r}_fmasked.nii -b fovea2.5dg+orig. -expr 'a*b' -prefix ${r}_fovea2.5_fmasked.nii
3dcalc -a ${r}_fmasked.nii -b periphery2.5dg+orig. -expr 'a*b' -prefix ${r}_periphery_fmasked.nii

echo 'V1, V2 & V3 fovea'
3dBrickStat -sum ${r}_fovea2.5_fmasked.nii;
echo 'V1, V2 & V3 periphery'
3dBrickStat -sum ${r}_periphery_fmasked.nii;
done

#make V2+3_fmasked.nii
3dmerge -gmax -prefix V2+3_fmasked.nii V2_fmasked.nii V3_fmasked.nii 

mv *_fmasked.nii ../../../ROImasks/


ROIs=( V1 V2 V3 )
# only V1, V2, V3 * fovea or periphery (no functional localizer thresholding!)
for r in "${ROIs[@]}"
do
3dcalc -a ${r}+orig. -b fovea2.5dg+orig. -expr 'a*b' -prefix ${r}_fovea2.5.nii
3dcalc -a ${r}+orig. -b periphery2.5dg+orig. -expr 'a*b' -prefix ${r}_periphery.nii

echo 'V1, V2 & V3 fovea'
3dBrickStat -sum ${r}_fovea2.5.nii;
echo 'V1, V2 & V3 periphery'
3dBrickStat -sum ${r}_periphery.nii;
done




####

ROIs=( fovea3.5dg para+fovea3.5dg )

for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix $r lh.${r}+orig. rh.${r}+orig.
done

ROIs=(fovea3.5dg)
for r in "${ROIs[@]}"
do
#seperate parafovea mask -> exclude the foveal mask part from parafoveal mask
3dcalc -a para+fovea3.5dg+orig. -b fovea3.5dg+orig. -expr 'a*b' -prefix tmp.overlap_${r}_peri
3dcalc -a para+fovea3.5dg+orig. -b tmp.overlap_${r}_peri+orig. -expr 'a-b' -prefix periphery-${r}

#seperate fovea mask -> exclude fixation part from foveal mask
#3dcalc -a fovea3.5dg+orig. -b fix0.61dg+orig. -expr 'a*b' -prefix tmp.overlap_${r}_fix
#3dcalc -a fovea3.5dg+orig. -b tmp.overlap_${r}_fix+orig. -expr 'a-b' -prefix ${r}-fix0.61dg
done

rm -f tmp.* 

ROIs=( V1 V2 V3 )
# threshold with functional contrast map
for r in "${ROIs[@]}"
do
3dcalc -a ${r}_fmasked.nii -b fovea3.5dg+orig. -expr 'a*b' -prefix ${r}_fovea3.5_fmasked.nii
#3dcalc -a ${r}_fmasked.nii -b fovea3.5dg-fix0.61dg+orig. -expr 'a*b' -prefix ${r}_fovea3.5-fix0.61_fmasked.nii
3dcalc -a ${r}_fmasked.nii -b periphery-fovea3.5dg+orig. -expr 'a*b' -prefix ${r}_periphery-fovea3.5_fmasked.nii
done


ROIs=( V1 V2 V3 )
# only V1, V2, V3 * fovea or periphery (no functional localizer thresholding!)
for r in "${ROIs[@]}"
do
3dcalc -a ${r}+orig. -b fovea3.5dg+orig. -expr 'a*b' -prefix ${r}_fovea3.5.nii
#3dcalc -a ${r}+orig. -b fovea3.5dg-fix0.61dg+orig. -expr 'a*b' -prefix ${r}_fovea3.5-fix0.61.nii
3dcalc -a ${r}+orig. -b periphery-fovea3.5dg+orig. -expr 'a*b' -prefix ${r}_periphery-fovea3.5.nii

done

ROIs=( V1 V2 V3 )
#fov+peri 
for r in "${ROIs[@]}"
do
3dmerge -gmax -prefix ${r}_fov+peri.nii ${r}_fovea3.5.nii ${r}_periphery-fovea3.5.nii 
done

mv *fovea3.5*.nii ../../../ROImasks
mv *fovea3.5*_fmasked.nii ../../../ROImasks/
mv ./*.nii ../../../ROImasks/

#####
    




 3dSurf2Vol -spec $surfdir/${surfID}_$1.spec \
  -surf_A $surfdir/$1.smoothwm.asc \
  -surf_B $surfdir/$1.pial.asc \
  -sv $svfile \
  -grid_parent $basegrid \
  -map_func max \
  -f_steps 10 \
  -f_p1_fr      -0.1  \
  -f_pn_fr      -0.1  \
  -f_index voxels  \
  -sdata_1D ${ROI_dset[$i]} \
  -prefix ${prefix[$i]}








