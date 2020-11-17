
MNI_DIR=/sas2/PECON/HJY/CrM/MNItemplate
tete=MNI152_T1_1mm_brain.nii
MNI2mm=MNI152_T1_2mm_brain.nii
standard=MNI152_T1_1mm_brain.nii
name=CIN@MNI

SNs=(02YYH 03LJA 04LSY 05YJG 06JDH 07KHY 08LKY 09JEH 10NJS 11KDB 12SCL 13NHJ 14JJY) 08LKY 
SN=01HJH

for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep
cp ../../mainPrep/CrM_CIN*.nii .
cp ../../mainPrep/T1_CIN.nii .

done

#s08 different structure 


#### Align T1 onto MNI space #### 
for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep

mv T1_CIN.nii 1_T1_CIN.nii 

#anatomical data 
3drefit -deoblique 1_T1_CIN.nii 

# change orientation of T1: from RAI to RPI (MNI space orientation)
3dresample -oreint RPI -prefix 2_rpi.T1_CIN.nii -inset 1_T1_CIN.nii

# Bias field correction
3dUnifize -input 2_rpi.T1_CIN.nii -prefix 3_unf.T1_CIN.nii -clfrac 0.5
3dSkullStrip -orig_vol -input 3_unf.T1_CIN.nii -prefix 4_SS.T1_CIN.nii

#register T1 to the template.
#Theoretically, dof 6 should works but dof 12 works better... 
flirt -in 4_SS.T1_CIN.nii -ref ${MNI_DIR}/${tete} -out T1_CIN@MNI2.nii -omat tmp.transf_HR2STD.mat -dof 12 -cost normmi


done

for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep

#volume registration using afni_proc.py(no alignment here)
name=CIN@MNI
dname=CrM_CIN
#CrM_vOnly
afni_proc.py -subj_id ${name} \
-dsets ${dname}1.nii ${dname}2.nii ${dname}3.nii ${dname}4.nii ${dname}5.nii \
 ${dname}6.nii ${dname}7.nii ${dname}8.nii ${dname}9.nii ${dname}10.nii \
-volreg_base_dset ${dname}5.nii['0'] \
-regress_censor_motion 0.5 \
-regress_censor_outliers 0.1 \
-blocks volreg mask regress
tcsh -xef proc.${name} |& tee output.proc.${name}
done


for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep


cd ${name}.results/

#copy T1 & template. we need these things later  
cp ../4_SS.T1_CIN.nii ../tmp.transf_*.mat . 
fslmaths ${MNI_DIR}/${tete} standard.nii

# change orientation from RAI to RPI (MNI space orientation)
run=(01 02 03 04 05 06 07 08 09 10)
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


done


for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep


cd ${name}.results/
rm -rf transf_*.mat
cp ../tmp.transf_*.mat . 
flirt -ref 4_SS.T1_CIN.nii -in refEPI.nii -out EPI@T1.nii -omat tmp.transf_EPI2HR.mat \
-cost normcorr -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear

# Read from specified file, or from standard input
infile="${1:-/sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep/CIN@MNI.results/tmp.transf_EPI2HR.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_EPI2HR.mat

# Read from specified file, or from standard input
infile="${1:-/sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep/CIN@MNI.results/tmp.transf_HR2STD.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_HR2STD.mat


###
convert_xfm -omat tmp.transf_EPI2STD.mat -concat transf_HR2STD.mat transf_EPI2HR.mat #concat two mats. why? ok I got it 

# Read from specified file, or from standard input
infile="${1:-/sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep/CIN@MNI.results/tmp.transf_EPI2STD.mat}"

while read line; do

    for number in $line; do
        printf "%f " "$number"
    done
    echo
done < $infile > transf_EPI2STD.mat

rm -f tmp.transf_*

#check first. Align EPI to MNI using transf_EPI2STD (use EPI2STD.mat here for the first guess)
flirt -ref ${MNI_DIR}/${MNI2mm} -in refEPI.nii -out refEPI@MNI.nii \
-cost corratio -applyxfm -init transf_EPI2STD.mat -interp trilinear

done


for SN in "${SNs[@]}"
do

cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep


cd ${name}.results/

#using ref.matrix, align the other images too (pb03.reg is the latest ones!)
for mm in "${run[@]}"
	do
	flirt -applyxfm -init transf_EPI2STD.mat -in pb02.${name}.r${mm}.RPI.nii \
	-ref ${MNI_DIR}/${MNI2mm} -out pb03.${name}.r${mm}.reg.nii -interp trilinear
done

#Align full mask@MNI
3drefit -deoblique full_mask.${name}+orig
3dresample -oreint RPI -prefix rpi.full_mask.${name}.nii -inset full_mask.${name}+orig

flirt -applyxfm -init transf_EPI2STD.mat -in rpi.full_mask.${name}.nii \
	-ref ${MNI_DIR}/${MNI2mm}.nii -out full_mask.${name}_reg -interp trilinear

run=(01 02 03 04 05 06 07 08 09 10)

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
	        -expr 'a+b' -prefix pb07.${name}.r${i}.sc_dt_hp_am  #for GLM
	done
	##Finished. Move the data -------------------------------------
	3dTcat -prefix ${SN}${name}_FNL.nii \
	pb07.${name}.r*.HEAD

	mv ${SN}${name}_FNL.nii ../../Searchlight_${SN}/
	
cp full_mask.${name}_reg.nii.gz ../../Searchlight_${SN}
	mv full_mask.${name}_reg.nii.gz ../../ROImasks
done



for SN in "${SNs[@]}"
do


cd /sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep/${name}.results
BH_DIR=/sas2/PECON/HJY/CrM/Exp2/${SN}/BH_data/onset

subj=CIN@MNI

3dcopy ../../Searchlight_${SN}/full_mask.${name}_reg.nii.gz ./full_mask.${name}_reg+tlrc.
#-------------- pb08 smoothing  ------------------# 

for i in "${run[@]}"
do
3dmerge -1blur_fwhm 5 -doall -prefix rm.pb08.${name}.r${i}.sc_dt_hp_am_blur pb07.${name}.r${i}.sc_dt_hp_am+tlrc.

# and set boundaries using anat mask 
3dcalc -a rm.pb08.${name}.r${i}.sc_dt_hp_am_blur+tlrc -b full_mask.${name}_reg+tlrc. \
-expr 'a*b' -prefix pb08.${name}.r${i}.sc_dt_hp_am_blur
done

rm -f rm.pb08*

#-------------- Regress ------------------# 
3dDeconvolve -input pb08.${name}.r*.sc_dt_hp_am_blur+tlrc.HEAD           \
    -censor censor_${name}_combined_2.1D                                    \
    -polort 3 							            \
    -num_stimts 8 			                                    \
    -stim_file 1 motion_demean.1D'[0]' -stim_base 1 -stim_label 1 roll      \
    -stim_file 2 motion_demean.1D'[1]' -stim_base 2 -stim_label 2 pitch     \
    -stim_file 3 motion_demean.1D'[2]' -stim_base 3 -stim_label 3 yaw       \
    -stim_file 4 motion_demean.1D'[3]' -stim_base 4 -stim_label 4 dS        \
    -stim_file 5 motion_demean.1D'[4]' -stim_base 5 -stim_label 5 dL     \
    -stim_file 6 motion_demean.1D'[5]' -stim_base 6 -stim_label 6 dP     \
    -stim_times 7 ${BH_DIR}/${SN}_AVonset.txt 'BLOCK(6,1)'           \
    -stim_label 7 multimodal                                                   \
    -stim_times 8 ${BH_DIR}/${SN}_Fixonset.txt 'BLOCK(8,1)'           \
    -stim_label 8 fix                                                \
    -local_times                                                            \
    -gltsym 'SYM: +multimodal -fix'                     \
    -glt_label 1 multimodal_vs_fix                      \
    -gltsym 'SYM: +fix -multimodal'                     \
    -glt_label 2 fix_vs_multimodal                      \
    -float                                                                  \
    -jobs 8                                                                 \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                 \
    -x1D_uncensored X.nocensore54.xmat.1D                                      \
    -bucket stats.${name}_multimodal -overwrite

done
done

#-------------------------------------------------

cd /sas2/PECON/HJY/CrM/Exp2/groupAnalysis

for SN in "${SNs[@]}"
do


MAIN_DIR=/sas2/PECON/HJY/CrM/Exp2/${SN}/Normalization/mainPrep/${name}.results


3dbucket -prefix subj_${SN}.beta+tlrc ${MAIN_DIR}/stats.${name}_multimodal+tlrc.HEAD'[multimodal#0_Coef, fix#0_Coef]'

done


3dttest++ -prefix all_t -brickwise -setA subj_*.beta+tlrc.HEAD -mask Mask_MNI+tlrc.


3drefit -sublabel 0 'Print_Coef' all_t+tlrc.
3drefit -sublabel 1 'Print_Tstat' all_t+tlrc.

3dttest++ -prefix all_t_clustsim  -Clustsim 12 \
-setA subj_*.beta+tlrc.HEAD'[multimodal#0_Coef]' \
-setB subj_*.beta+tlrc.HEAD'[fix#0_Coef]'\
-mask Mask_MNI+tlrc.


	cp pb07.${name}.r*.nii	/sas2/PECON/HJY/CrM/Exp2/groupAnalysis

















