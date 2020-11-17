#!/bin/bash
#for Crossmodal experiment 2: 
#
#session1: CrM_VA (114TR), visLoc, & audLoc
#session2: CrM_CIN (148TR)
#
#CrM_CIN EPI image closest to T1 will be the norm!
#This time, run this script after acquiring all session data
#all data including localizer ones should be in the 'mainPrep' folder
#only raw data and final data keep .nii 

surfID=190906KDB
SN=11KDB
session= #CrM_VA or CRM_CIN

SURF_DIR=/sas2/PECON/freesurfer/subjects
CRM_DIR=/sas2/PECON/HJY/CrM/Exp2/${SN}/mainPrep

cd ${CRM_DIR}
##------------------- CrM_CIN (session 2) ALIGNMENT ------------------## 
seNm=CIN #session name

#Session2 EPI data will be the norm! 
3dUnifize -input T1_${seNm}_SS+orig. -prefix T1_${seNm}_SS_uni -clfrac 0.5 

#align Inplane1 to EPI first (Make Inplane1@CrM)
@Align_Centers -no_cp -base CrM_${seNm}4.nii -dset inplane_${seNm}_SS+orig.
align_epi_anat.py -dset1 CrM_${seNm}4.nii[2] -dset2 inplane_${seNm}_SS+orig. -dset2to1 \
-cost lpa -deoblique off -feature_size 0.5  -ginormous_move 
3drename inplane_${seNm}_SS_al+orig. Inplane_${seNm}@CrMvr2   #CrM version 2

#align T1 to Inplane@CrMvr2 (Make AVol@CrM)
@Align_Centers -no_cp -base Inplane_${seNm}@CrMvr2+orig. -dset T1_${seNm}_SS_uni+orig.
align_epi_anat.py -dset1 Inplane_${seNm}@CrMvr2+orig. -dset2 T1_${seNm}_SS_uni+orig. -dset2to1 \
-cost mi -deoblique off -feature_size 0.5  -ginormous_move 
3drename T1_${seNm}_SS_uni_al+orig. AVol@CrMvr2

#align surface volume to AVol@CrMvr2 (Make SVol@CrMvr2)
#choose options! 
@SUMA_AlignToExperiment -exp_anat AVol@CrMvr2+orig \
	-surf_anat $SURF_DIR/${surfID}/SUMA/${surfID}_SurfVol_SS+orig. \
	-strip_skull neither \
	-align_centers -al_opt 'cost mi' \
	-prefix SVol@CrMvr2 #sometimes it's better to leave -wd option out  # -wd -al_opt 'cost lpa+' \

#motion correction & alignment
afni_proc.py -subj_id CrM_${seNm} \
-dsets CrM_${seNm}1.nii CrM_${seNm}2.nii CrM_${seNm}3.nii CrM_${seNm}4.nii \
	CrM_${seNm}5.nii CrM_${seNm}6.nii CrM_${seNm}7.nii CrM_${seNm}8.nii \
	CrM_${seNm}9.nii CrM_${seNm}10.nii \
-volreg_base_dset CrM_${seNm}4.nii['2'] \
-copy_anat AVol@CrMvr2+orig. \
-volreg_align_e2a \
-regress_censor_motion 0.5 \
-regress_censor_outliers 0.1 \
-blocks align volreg mask regress
tcsh -xf proc.CrM_${seNm} 


##------------------- CrM_VAsep (session 3 which is an additional one) ALIGNMENT ------------------## 
seNm3=aOnly #session name

## CrM_VA
#align Inplane2 images to AVol@CrM+orig. or Inplane2@CrM+orig.

#Skull stripping 
3dcopy T1_${seNm3}.nii T1_${seNm3}
3dSkullStrip -input inplane_${seNm3}.nii -prefix inplane_${seNm3}_SS
3dUnifize -input AVol@CrMvr2+orig. -prefix AVol@CrMvr2_uni -clfrac 0.5 

#align session3 Inplane to AVol@CrMvr2 
@Align_Centers -no_cp -base AVol@CrMvr2_uni+orig. -dset inplane_${seNm3}_SS+orig.
align_epi_anat.py -dset1 AVol@CrMvr2_uni+orig. -dset2 inplane_${seNm3}_SS+orig. -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move
3drename inplane_${seNm3}_SS_al Inplane_${seNm3}@CrMvr2

#align session3 T1 to session1 AVol@CrMvr2 (if you have T1...)
3dSkullStrip -input T1_${seNm3}+orig. -prefix T1_${seNm3}_SS
3dUnifize -input T1_${seNm3}_SS+orig. -prefix T1_${seNm3}_SS_uni -clfrac 0.5 
@Align_Centers -no_cp -base AVol@CrMvr2+orig. -dset T1_${seNm3}_SS_uni+orig.
align_epi_anat.py -dset1 AVol@CrMvr2+orig. -dset2 T1_${seNm3}_SS_uni+orig. -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move
3drename T1_${seNm3}_SS_uni_al+orig. AVol_${seNm3}@CrMvr2

run=( 1 2 3 4 5 6 7 8 9 10 )  #CrM_aOnly
folder=(CrM_vOnly)  #CrM_aOnly  VAsep_aOnly

for ww in "${folder[@]}"
do

echo $ww
echo "${run[*]}"



for mm in "${run[@]}"
do
@Align_Centers -no_cp -base Inplane_${seNm3}@CrMvr2+orig. -dset CrM_${ww}${mm}.nii
done

#motion correction & alignment
afni_proc.py -subj_id CrM_aOnly \
-dsets CrM_${ww}1.nii CrM_${ww}2.nii CrM_${ww}3.nii CrM_${ww}4.nii \
	CrM_${ww}5.nii CrM_${ww}6.nii CrM_${ww}7.nii CrM_${ww}8.nii \
	CrM_${ww}9.nii CrM_${ww}10.nii \
-copy_anat Inplane_${seNm3}@CrMvr2+orig. \
-volreg_base_dset CrM_${ww}1.nii['2'] \
-volreg_align_e2a \
-align_opts_aea -cost lpc+ZZ -giant_move \
-regress_censor_motion 0.5 \
-regress_censor_outliers 0.1 \
-blocks align volreg mask regress
tcsh -xf proc.CrM_aOnly

done
run=( 01 02 03 04 05 06 07 08 09 10 )  #CrM_aOnly
##------------------- CrM_loc (session ?) ALIGNMENT ------------------## 
#use Inplane_${seNm}@CrMvr2 here
seNm3=audLoc
v_seNm=$seNm2 #same session with visLoc
a_seNm=$seNm2 #same session with audLoc

#In case you have different session for loc: Inplane to AVol@CrMvr2 
3dSkullStrip -input T1_CrM_${seNm3}.nii -prefix T1_${seNm3}_SS
@Align_Centers -no_cp -base AVol@CrMvr2+orig. -dset T1_${seNm3}_SS+orig.
align_epi_anat.py -dset1 AVol@CrMvr2+orig. -dset2 T1_${seNm3}_SS+orig. -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move
3drename T1_${seNm3}_SS_al AVol_${seNm3}@CrMvr2

3dSkullStrip -input inplane_CrM_${seNm3}.nii -prefix inplane_${seNm3}_SS
@Align_Centers -no_cp -base AVol@CrMvr2+orig. -dset inplane_${seNm3}_SS+orig.
align_epi_anat.py -dset1 AVol@CrMvr2+orig. -dset2 inplane_${seNm3}_SS+orig. -dset2to1 \
 -cost mi -deoblique off -feature_size 0.5 -ginormous_move
3drename inplane_${seNm3}_SS_al Inplane_${seNm3}@CrMvr2

mv *Loc*.nii ../localizerAnalysis
cp AVol_VA@CrMvr2+orig.* ../localizerAnalysis
cp Inplane*@CrMvr2+orig.* ../localizerAnalysis
cp AVol@CrMvr2+orig.* ../localizerAnalysis
cp SVol@CrMvr2+orig.* ../localizerAnalysis


##-------------------PREPROCESSING------------------## 
##-------------------PREPROCESSING------------------## 
##-------------------PREPROCESSING------------------## 
run=( 01 02 03 04 05 06 07 08 09 10 )
folder=(CrM_VAsep_vOnly)  #CrM_VAsep_vOnly CrM_VAsep_aOnly CrM_CIN

for ww in "${folder[@]}"
do

if [ "$ww" == "CrM_VA" ]; then
run=( 01 02 03 04 05 06 07 08 09 10 )
elif [ "$ww" == "CrM_CIN" ]; then
run=( 01 02 03 04 05 06 07 08 09 10 )
fi
echo $ww
echo "${run[*]}"

cd ${CRM_DIR}/${ww}.results

#gedit /home/jiyeongha/Script/phaseCorr.sh &

	##------------------------02 Scale & remove skull in EPI -----------

	for r in "${run[@]}"
	do
	3dClipLevel pb01.${ww}.r${r}.volreg+orig.HEAD >> clip.txt
	3dTstat -mean -prefix r.${r}.base pb01.${ww}.r${r}.volreg+orig.HEAD'[0..$]'
	done
	more clip.txt # check the smallest clip value across all runs
	clip=$(cut -f1 -d"," clip.txt | sort -n | head -1) # assign the clip value

	for r in "${run[@]}"
	do
	3dcalc -a pb01.${ww}.r${r}.volreg+orig. -b r.${r}.base+orig. \
	       -expr "(100 * a/b) * step(b-$clip)" -prefix pb02.${ww}.r${r}.scaled #remove the space and scale 
	done
	# (tip) when calling variables(like $clip) within 3dcalc -expr, be sure to use "", and not ''

	##------------------------ 03 Detrending & 04 highpass filter & 05 add mean (100)-----
	for i in "${run[@]}"
	do

	3dDetrend -polort 1 -prefix pb03.${ww}.r${i}.sc_dt pb02.${ww}.r${i}.scaled+orig.
	3dBandpass -prefix pb04.${ww}.r${i}.sc_dt_hp 0.01 99999 pb03.${ww}.r${i}.sc_dt+orig
	#add mean of scaled data (pb02.~)
	3dTstat -mean -prefix r.${i}.sc_base pb02.${ww}.r${i}.scaled+orig.HEAD'[0..$]'
	3dcalc  -a pb04.${ww}.r${i}.sc_dt_hp+orig.HEAD -b r.${i}.sc_base+orig.HEAD \
	        -expr 'a+b' -prefix pb05.${ww}.r${i}.sc_dt_hp_am
	done
	##Finished. Let's combine all data-------------------------------------
	3dTcat -prefix ${SN}${ww}_FNL.nii \
	pb05.${ww}.*.HEAD 
	mv ${SN}${ww}_FNL.nii ../../Decoding_${SN}/
done


SNs=( 01HJH 02YYH 03LJA 04LSY 05YJG 06JDH 07KHY 10NJS 11KDB )

for SN in "${SNs[@]}"
do
cp -r /sas2/PECON/HJY/CrM/Exp2/newAlignment_EPI@AVol/${SN}/Decoding_${SN}_newAlign /home/jiyeongha/Decoding_Exp2_newAlign_EPI@AVol/
done
















