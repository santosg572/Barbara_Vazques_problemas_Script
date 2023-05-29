santosg@geminis:~/Barbara_Vazques_problemas_Script$ cat scripts_bids_GLF2.sh 

#Archivos de lap a clúster
#rsync -e /home/lau/Documentos/longitudinal/GL_12M navalj@penfield:/misc/cannabis/alcauter/LauraMonica/longitudinal/


#!/bin/bash

sourceDir=/misc/geminis/navalj/GLF/

############# CHECK for the parameters, files and PATH specified in this SECTION ###################

TR=1;
sub=*
#sub=sub-24MD1s94

# In Hz
hpf=0.01 
lpf=0.08

# Reference Volume for motion correction
refvol=0;

# ATLAS FILES  ### Check for the correct altasdir and atlas files !!! and their PATH
atlasdir=/misc/geminis/navalj/prueba/atlas/TohokuUniv/correctOrientation; 
atlasbrain=${atlasdir}/brain.nii.gz;
atlasbrain5mm=${atlasdir}/brain_5mm.nii.gz;
atlasbrain5mmMask=${atlasdir}/brain_5mm_mask.nii.gz;

#  Voxel size (~x10)  ###  CORRECT FOR YOUR DATA !!!. All dimensions must be > 1, but keep the scale
mkdir -p ${sourceDir}/derivatives/ppBOLD
ls -d ${sourceDir}/${sub} | parallel cp -r {} ${sourceDir}/derivatives/ppBOLD

ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*gz | parallel fslchpixdim  {} 4.68 4.68 12
ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/*gz | parallel fslchpixdim  {} 1.17 1.17 12

############################  PREPROCESSING ########################################################

############# Slice timing and motion correction

ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold.nii.gz | cut -d . -f 1 | parallel slicetimer -i {} -o {}_pp --odd 
ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp.nii.gz | parallel mcflirt -in {} -refvol $refvol -plots

mkdir -p ${sourceDir}/derivatives/QC_pp
slicesdir $(ls ${sourceDir}/derivatives/ppBOLD/*/func/sub-*_task-rest_bold_pp.nii.gz)
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_func

############# CORREGISTROS    

### BET Estructural
ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w.nii.gz | cut -d . -f 1 | parallel N4BiasFieldCorrection -i {}.nii.gz -o {}_pp.nii.gz -t [0.3,0.01,100];
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w_p.nii.gz | cut -d . -f 1 | parallel DenoiseImage -i {}.nii.gz -o {}p.nii.gz;
ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w.nii.gz | cut -d . -f 1 | parallel fslcpgeom {} {}_pp;
ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w_pp.nii.gz | parallel fslswapdim {} x y z {}

for i in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w_pp.nii.gz);do
	echo "Processing " $(basename $i .nii.gz | cut -d _ -f 1); 
	fslchpixdim  $i 1.17 1.17 6
	j=$(echo $i | cut -d . -f 1);
	center=$(cluster -i $i -t 10 | sed -n '2{p;q}' | awk '{print int($7)" "int($8)" "int($9)}');
	bet ${i} ${j}_brain.nii.gz -r 75 -c $center -f 0.3 -g 0.25;
	bet ${j}_brain.nii.gz ${j}_brain.nii.gz -r 75 -c $center -f 0.25 -g -0.25;

	fslchpixdim  $i 1.17 1.17 12
	fslchpixdim ${j}_brain.nii.gz  1.17 1.17 12

done

slicesdir -o $(ls -r ${sourceDir}/derivatives/ppBOLD/*/anat/sub-*_T2w_pp*.nii.gz)
rm -rf slicesdir_BET; mv slicesdir/ ${sourceDir}/derivatives/QC_pp/slicesdir_BET

### T2 brain to ATLAS

# register T2brain to atlas
ls  ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w_pp_brain.nii.gz | cut -d . -f 1 | parallel flirt -in {} -ref $atlasbrain -out {}_2TohokuA -omat {}_2TohokuA.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear;

# Checar corregistros al atlas
slicesdir -p $atlasbrain $(ls ${sourceDir}/derivatives/ppBOLD/*/anat/sub-*_T2w_pp_brain_2TohokuA.nii.gz)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_T2toAtlas; mv slicesdir/ ${sourceDir}/derivatives/QC_pp/slicesdir_T2toAtlas

### rsfMRI to T2

# reorient rsfMRIs
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp.nii.gz | cut -d . -f 1 | parallel fslswapdim {} x y z {}
# ExampleFunc
ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp.nii.gz | cut -d . -f 1 | parallel fslroi {} {}_examplefunc $refvol 1;
# rsfMRI to T2
for i in $(ls ${sourceDir}/derivatives/ppBOLD/${sub}/anat/sub-*_T2w_pp.nii.gz);do
	id=$(basename $i | cut -d _ -f 1)
	echo "Processing "$i;
	for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/func/${id}*_examplefunc.nii.gz | cut -d . -f 1);do 
	echo $(basename $j);
flirt -in $j -ref $i -out ${j}_2T2 -omat ${j}_2T2.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -2D -dof 7 -interp trilinear;
#-searchrx -90 90 -searchry -90 90 -searchrz -90 90 -2D 
	done
done

rm ${sourceDir}/derivatives/listcheck_exfunc2T2.list
for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do id=$(basename $i);echo ${sourceDir}/derivatives/ppBOLD/${id}/func/${id}_task-rest_bold_pp_examplefunc_2T2.nii.gz >> ${sourceDir}/derivatives/listcheck_exfunc2T2.list;echo ${sourceDir}/derivatives/ppBOLD/${id}/anat/${id}_T2w.nii.gz >> ${sourceDir}/derivatives/listcheck_exfunc2T2.list;done
slicesdir -o $(cat ${sourceDir}/derivatives/listcheck_exfunc2T2.list)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_EF2T2w
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_EF2T2w

### rsfMRI to ATLAS

# Combine transformations and apply

for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do echo $i;
	id=$(basename $i);
	for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/func/${id}*_bold.nii.gz | cut -d . -f 1);do 
	#combine
	convert_xfm -omat ${j}_pp_examplefunc_2atlasbrain.mat ${sourceDir}/derivatives/ppBOLD/${id}/anat/${id}_T2w_pp_brain_2TohokuA.mat -concat  ${j}_pp_examplefunc_2T2.mat;
	# Apply
	flirt -in ${j}_pp_mcf.nii.gz -applyxfm -init ${j}_pp_examplefunc_2atlasbrain.mat -out ${j}_pp_mcf_2atlasbrain5mm -paddingsize 0.0 -interp trilinear -ref $atlasbrain5mm;
	done	
done

slicesdir -p $atlasbrain5mm $(ls ${sourceDir}/derivatives/ppBOLD/*/func/*_pp_mcf_2atlasbrain5mm.nii.gz)
rm -rf ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm
mv slicesdir ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm

############# aCompCor and Band Passing 

for i in $(ls -d ${sourceDir}/derivatives/ppBOLD/${sub});do echo $i;
	id=$(basename $i);
	for j in $(ls ${sourceDir}/derivatives/ppBOLD/${id}/func/${id}*_bold_pp_mcf_2atlasbrain5mm.nii.gz | cut -d . -f 1);do

	liminf=$(fslstats $j -r | awk '{print $1}'); liminf2=$(echo "scale=2; 2*${liminf}" | bc);
	${FSLDIR}/bin/fslmaths $j -Tmean -thr $liminf2 -bin -mas ${atlasdir}/CSFWMmask_5mm.nii.gz ${j}_CSFWMmask_5mm.nii.gz
	${FSLDIR}/bin/fslmeants -i $j -o ${j}_aCompCor.txt -m ${j}_CSFWMmask_5mm.nii.gz  --eig --order=5;

	rm ${j}_aCCtf.nii.gz
	3dTproject -input ${j}.nii.gz -prefix ${j}_aCCtf.nii.gz -polort 0 -ort ${j}_aCompCor.txt -passband $hpf $lpf -TR ${TR} -mask $atlasbrain5mmMask

	done
done

### Smoothing FWHM 10mm sigma=4.25
ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*bold_pp_mcf_2atlasbrain5mm_aCCtf.nii.gz | cut -d . -f 1 | parallel fslmaths {} -s 4.25 -mas $atlasbrain5mmMask {}_sm425.nii.gz

############# MENSAJE FINAL #############

#firefox --new-window ${sourceDir}/derivatives/QC_pp/slicesdir_func/index.html;
#firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_BET/index.html;
#firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_T2toAtlas/index.html;
#firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_EF2T2w/index.html;
#firefox --new-tab ${sourceDir}/derivatives/QC_pp/slicesdir_ppAtlas5mm/index.html;
#echo "        DONE"
#echo ""


############# Optional but RECOMMENDED to save space:
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp_mcf.nii.gz
# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/sub-*_task-rest_bold_pp_mcf_2atlasbrain5mm.nii.gz



################## TimeFiltering and smoothing Alternatives (using FSL only) #########################

### Smoothing: DenoiseImage (takes about 2-3 hours per 4D series)
#ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*bold_pp_mcf_2atlasbrain5mm_aCCtf.nii.gz | cut -d . -f 1 | parallel DenoiseImage -i {} -x $atlasbrain5mmMask -o {}_smDI.nii.gz

### aCompCor and time filtering with FSL

# fsl_glm
# ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*bold_pp_mcf_2atlasbrain5mm.nii.gz | cut -d . -f 1 | parallel fsl_glm -i {} -m $atlasdirbrain5mmMask -d {}_aCompCor.txt --demean --out_res={}_aCC.nii.gz;

# TR=1;
# hpf=0.01;
# lpf=0.1;
# hp_sigma=`echo "scale=2 ;(1/${hpf})/2.35/${TR}" | bc`; # In volumes for fslmaths
# lp_sigma=`echo "scale=2 ;(1/${lpf})/2.35/${TR}" | bc`; # In volumes for fslmaths

# ls ${sourceDir}/derivatives/ppBOLD/${sub}/func/*bold_pp_mcf_2atlasbrain5mm_aCC.nii.gz | cut -d . -f 1 | parallel fslmaths {} -bptf  $hp_sigma $lp_sigma -mas $atlasbrain5mmMask {}tf.nii.gz

# rm ${sourceDir}/derivatives/ppBOLD/${sub}/func/*bold_pp_mcf_2atlasbrain5mm_aCC.nii.gz






santosg@geminis:~/Barbara_Vazques_problemas_Script$ 
