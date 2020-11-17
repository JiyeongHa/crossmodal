
for ((k=1;k<=13;k++));
do

SNs=(01LES 02LJA 03HJJ 04HJH 05YYH 06YJG 07CES 08LKY 09JHY 10LIY 11KDB 12NJS 13NHJ 14JJY)
retIDs=(133 122 136 123 132 64 128 143 144 126 180 178 188 181)
ver=ret

retID=${retIDs[$k]}
SN=${SNs[$k]}


mkdir -p /sas2/PECON/HJY/CrM/Exp1/groupAnalysis/eccenBias/${SN}_fromRet/mainPrep
mkdir -p /sas2/PECON/HJY/CrM/Exp1/groupAnalysis/eccenBias/${SN}_fromRet/ret@MNI_${SN}
MAIN_DIR=/sas2/PECON/HJY/CrM/Exp1/groupAnalysis/eccenBias/${SN}_fromRet/mainPrep
RET_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/Img_data/
MNI_DIR=/sas2/PECON/HJY/CrM/MNItemplate
tete=MNI152_T1_1mm_brain.nii
MNI2mm=MNI152_T1_2mm_brain.nii
cd ${MAIN_DIR}

cp ${RET_DIR}/ring1.nii .
cp ${RET_DIR}/ring2.nii .
cp ${RET_DIR}/T1.nii 1_T1_${ver}.nii 


#### Align T1 onto MNI space #### 

#anatomical data 
3drefit -deoblique 1_T1_${ver}.nii 

# change orientation of T1: from RAI to RPI (MNI space orientation)
3dresample -oreint RPI -prefix 2_rpi.T1_${ver}.nii -inset 1_T1_${ver}.nii

# Bias field correction
3dUnifize -input 2_rpi.T1_${ver}.nii -prefix 3_unf.T1_${ver}.nii -clfrac 0.5
3dSkullStrip -orig_vol -input 3_unf.T1_${ver}.nii -prefix 4_T1_SS.${ver}.nii

#register T1 to the template.
#Theoretically, dof 6 should works but dof 12 works better... 
flirt -in 4_T1_SS.${ver}.nii -ref ${MNI_DIR}/${tete} -out T1@MNI.nii -omat tmp.transf_HR2STD.mat -dof 12 -cost normmi 

#volume registration using afni_proc.py(no alignment here)
name=${ver}@MNI
dname=ring
#CrM_vOnly
afni_proc.py -subj_id ${name} \
-dsets ${dname}1.nii ${dname}2.nii \
-volreg_base_dset ${dname}1.nii['0'] \
-regress_censor_motion 0.5 \
-regress_censor_outliers 0.1 \
-blocks volreg mask regress
tcsh -xef proc.${name} |& tee output.proc.${name}

cd ${name}.results/

#copy T1 & template. we need these things later  
cp ../4_*SS*.nii ../tmp.transf_*.mat . 
fslmaths ${MNI_DIR}/${tete} standard.nii

# change orientation from RAI to RPI (MNI space orientation)
run=(01 02)
for mm in "${run[@]}"
do
fslmaths -odt float #(do we need this?)
3drefit -deoblique pb01.${name}.r${mm}.volreg+orig.
3dresample -oreint RPI -prefix pb02.${name}.r${mm}.RPI.nii -inset pb01.${name}.r${mm}.volreg+orig.
done

#Make ref. EPI: 1TR copy from vol-regged EPI data
#It'll be used as an reference image, as in 
#ref@MNI and then use the output matrix to the other EPI images. 
fslroi pb02.${name}.r${mm}.RPI.nii refEPI.nii 0 -1 0 -1 0 -1 0 1
flirt -ref 4_T1_SS.${ver}.nii -in refEPI.nii -out EPI@T1.nii -omat tmp.transf_EPI2HR.mat \
-cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear

# Read from specified file, or from standard input
infile="${1:-tmp.transf_EPI2HR.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_EPI2HR.mat

# Read from specified file, or from standard input
infile="${1:-tmp.transf_HR2STD.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_HR2STD.mat


###
convert_xfm -omat tmp.transf_EPI2STD.mat -concat transf_HR2STD.mat transf_EPI2HR.mat #concat two mats. why? ok I got it 

# Read from specified file, or from standard input
infile="${1:-tmp.transf_EPI2STD.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_EPI2STD.mat

#rm -f tmp.transf_*

#check first. Align EPI to MNI using transf_EPI2STD (use EPI2STD.mat here for the first guess)
flirt -ref ${MNI_DIR}/${MNI2mm} -in refEPI.nii -out refEPI@MNI.nii \
-cost corratio -applyxfm -init transf_EPI2STD.mat -interp trilinear

#using ref.matrix, align the other images too (pb03.reg is the latest ones!)
for mm in "${run[@]}"
	do
	flirt -applyxfm -init transf_EPI2STD.mat -in pb02.${name}.r${mm}.RPI.nii \
	-ref ${MNI_DIR}/${MNI2mm} -out pb03.${name}.r${mm}.reg.nii -interp trilinear
done


	##------------------------04 Scale & remove skull in EPI -----------

	for r in "${run[@]}"
	do
	3dClipLevel pb03.${name}.r${r}.reg.nii >> clip.txt
	3dTstat -mean -prefix r.${r}.base.nii pb03.${name}.r${r}.reg.nii'[0..$]'
	done
	more clip.txt # check the smallest clip value across all runs
	clip=$(cut -f1 -d"," clip.txt | sort -n | head -1) # assign the clip value

	for r in "${run[@]}"
	do
	3dcalc -a pb03.${name}.r${r}.reg.nii -b r.${r}.base.nii \
	       -expr "(100 * a/b) * step(b-$clip)" -prefix pb04.${name}.r${r}.scaled.nii #remove the space and scale 
	done
	# (tip) when calling variables(like $clip) within 3dcalc -expr, be sure to use "", and not ''

	##------------------------ 05 Detrending & 06 highpass filter & 07 add mean (100)-----
	for i in "${run[@]}"
	do

	3dDetrend -polort 1 -prefix pb05.${name}.r${i}.sc_dt.nii pb04.${name}.r${i}.scaled.nii
	3dBandpass -prefix pb06.${name}.r${i}.sc_dt_hp.nii 0.01 99999 pb05.${name}.r${i}.sc_dt.nii
	#add mean of scaled data (pb02.~)
	3dTstat -mean -prefix r.${i}.sc_base.nii pb04.${name}.r${i}.scaled.nii'[0..$]'
	3dcalc  -a pb06.${name}.r${i}.sc_dt_hp.nii -b r.${i}.sc_base.nii \
	        -expr 'a+b' -prefix pb07.${name}.r${i}.sc_dt_hp_am.nii
	done
	##Finished. Let's combine all data-------------------------------------
	3dTcat -prefix ${SN}ring@MNI_FNL.nii \
	pb07.${name}.*.nii
	mv ${SN}ring@MNI_FNL.nii ../../ret@MNI_${SN}/





SNs=(01LES 02LJA 03HJJ 04HJH 05YYH 06YJG 07CES 08LKY 09JHY 10LIY 11KDB 12NJS 13NHJ 14JJY)
retIDs=(133 122 136 123 132 64 128 143 144 126 180 178 188 181)
ver=ret
subj=ret@MNI
for ((k=0;k<=13;k++));
do

retID=${retIDs[$k]}
SN=${SNs[$k]}
MAIN_DIR=/sas2/PECON/HJY/CrM/Exp1/groupAnalysis/eccenBias/${SN}_fromRet/mainPrep/ret@MNI.results
RET_DIR=/group_hpc/WMShimLab/Retinotopy/${retID}Retinotopy/Img_data/
MNI_DIR=/sas2/PECON/HJY/CrM/MNItemplate
STAT_DIR=/sas2/PECON/HJY/CrM/Exp1/groupAnalysis/eccenBias/GLM
tete=MNI152_T1_1mm_brain.nii
MNI2mm=MNI152_T1_2mm_brain.nii
cd ${MAIN_DIR}


for i in "${run[@]}"
do
3dmerge -1blur_fwhm 5 -doall -prefix rm.pb08.${subj}.r${i}.sc_dt_hp_am_blur    \
pb07.${subj}.r${i}.sc_dt_hp_am.nii

# and set boundaries using anat mask 
3dcalc -a rm.pb08.${subj}.r${i}.sc_dt_hp_am_blur+tlrc -b /sas2/PECON/HJY/CrM/Exp1/groupAnalysis/Mask_MNI+tlrc. \
-expr 'a*b' -prefix pb08.${subj}.r${i}.sc_dt_hp_am_blur+tlrc.
done

rm -f rm.pb08*

#-------------- Regress ------------------# 
3dDeconvolve -input pb08.${subj}.r*.sc_dt_hp_am_blur+tlrc.HEAD            \
    -censor censor_${subj}_combined_2.1D                                    \
    -polort 3 							            \
    -num_stimts 8 			                                    \
    -stim_file 1 motion_demean.1D'[0]' -stim_base 1 -stim_label 1 roll      \
    -stim_file 2 motion_demean.1D'[1]' -stim_base 2 -stim_label 2 pitch     \
    -stim_file 3 motion_demean.1D'[2]' -stim_base 3 -stim_label 3 yaw       \
    -stim_file 4 motion_demean.1D'[3]' -stim_base 4 -stim_label 4 dS        \
    -stim_file 5 motion_demean.1D'[4]' -stim_base 5 -stim_label 5 dL     \
    -stim_file 6 motion_demean.1D'[5]' -stim_base 6 -stim_label 6 dP     \
    -stim_times 7 /home/jiyeongha/Script/center_onset.txt 'BLOCK(2,1)'           \
    -stim_label 7 center                                                  \
    -stim_times 8 /home/jiyeongha/Script/periphery_onset.txt 'BLOCK(2,1)'           \
    -stim_label 8 periphery                                                \
    -local_times                                                            \
    -gltsym 'SYM: +center -periphery'                     \
    -glt_label 1 center_vs_periphery                      \
    -gltsym 'SYM: +periphery -center'                     \
    -glt_label 2 periphery_vs_center                      \
    -float                                                                  \
    -jobs 8                                                                 \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                 \
    -x1D_uncensored X.nocensore54.xmat.1D                                      \
    -bucket stats.${subj}_eccenBias

3dbucket -prefix ${STAT_DIR}/subj_${SN}.beta+tlrc stats.${subj}_eccenBias+tlrc.HEAD'[center#0_Coef, periphery#0_Coef]'
done



##ttest


3dttest++ -prefix all_t_pairwise -brickwise -setA subj_*.beta+tlrc.HEAD'[center#0_Coef]' \
-setB subj_*.beta+tlrc.HEAD'[periphery#0_Coef]' -mask /sas2/PECON/HJY/CrM/Exp1/groupAnalysis/Mask_MNI+tlrc. 

3drefit -sublabel 0 'Print_Coef' all_t+tlrc.
3drefit -sublabel 1 'Print_Tstat' all_t+tlrc.

3dttest++ -prefix all_t_clustsim  -Clustsim 12 \
-setA subj_*.beta+tlrc.HEAD'[Aud#0_Coef]' \
-setB subj_*.beta+tlrc.HEAD'[fix#0_Coef]'\
-mask ../Mask_MNI+tlrc.



3dttest++ -prefix all_t_clustsim.nii  -Clustsim 12 \
-setA subj_*.beta+tlrc.HEAD'[center#0_Coef]' \
-setB subj_*.beta+tlrc.HEAD'[periphery#0_Coef]'\
-mask /sas2/PECON/HJY/CrM/Exp1/groupAnalysis/Mask_MNI+tlrc. 

3dttest++ -prefix sl_one_t_clustsim  -Clustsim 12 \
-setA *_sl-0.5+tlrc.HEAD \
-mask ../../Mask_MNI+tlrc.

for SN in "${SNs[@]}"
do
#3dcopy /sas2/PECON/HJY/CrM/Exp1/${SN}/Normalization/Searchlight_${SN}/${SN}_sl.nii.gz ./${SN}_sl
3dcalc -a ./${SN}_sl+tlrc. -expr 'a-0.5' -prefix ./${SN}_sl-0.5
done






