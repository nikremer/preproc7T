%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Diffusion MRI preprocessing pipeline using MRtrix3 + FSL + some Matlab
%
%
% updated: 14 dec 2018 
% by: L. Liebrand
%
% (Incomplete list of) dependencies: Jimmy Shen toolbox, FSL 5.0.10, 
% MRtrix 3, bash scripts: dti_topup_eddy.sh, total_readout_time_dcm, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% This pipeline performs standard diffusion MRI preprocessing. It expects
% raw diffusion scan (nifti) + .bvec and .bval files, raw topup scan
% (nifti), 2 raw dicom files (data + topup scan) for reading topup variables
% dti_raw.nii, dti_raw.bvec, dti_raw.bval, dti_topup_raw.nii


%% Files
% Unzip dti_raw and dti_topup_raw %barbara
!gzip -d dti_raw.nii
!gzip -d dti_topup_raw.nii

% raw data filename
dti='dti_raw.nii';       % make sure it is not dti_raw.nii.gz                       
[~,dti,~]=fileparts(dti);
% raw topup filename
dti_tp='dti_topup_raw.nii';                      
[~,dti_tp,~]=fileparts(dti_tp);

% pipeline requires one dicom file of diffusion + topup series for readout
% time calculation
%dcm='/data/lcliebrand/scan_data/20190129_3T_hires/MR/00501_dwi56_1.28mm_b800/11653.dcm';
%dcm_tp='/data/lcliebrand/scan_data/20190129_3T_hires/MR/00601_dwi_1.28mm_b0_TOPUP/00299.dcm';

%% First preprocessing steps (gibbs artifacts, denoising etc)

% make directory for processing
if ~exist('dti7t','dir'); mkdir('dti7t'); end   

%ls('/usr/local/mrtrix3/bin')
% perform mrtrix denoising + de-gi bbs-ringing on raw data (mrtrix or script by dr. M.W.A. Caan)
unix(['/usr/local/mrtrix3/bin/dwidenoise -nthreads 14 -extent 3,3,3 -noise dti7t/dti_noisemap.nii ',dti,'.nii dti7t/dti_denoised.nii'])
unix(['/usr/local/mrtrix3/bin/mrcalc ',dti,'.nii dti7t/dti_denoised.nii -subtract dti7t/dti_residuals.nii'])
unix('/usr/local/mrtrix3/bin/mrdegibbs -nthreads 14 -axes 0,1 dti7t/dti_denoised.nii dti7t/dti_unring.nii')

%fixB0RingingArtifact('dti7t/dti_denoised.nii','dti_raw','dti7t/dti_unring.nii'); 

% perform mrtrix denoising + de-gibbs-ringing on raw topup data
% NB: IMPOSSIBLE if topup data is only 1 volume

[~,cmdout]=(unix('fslinfo dti_topup_raw.nii | grep ''dim4'' | awk ''NR==1{print$2}'' | tr -d  ''\n''    '));
if cmdout~='1'
    unix(['/usr/local/mrtrix3/bin/dwidenoise -nthreads 14 -noise dti7t/dti_topup_noisemap.nii ',dti_tp,'.nii dti7t/dti_topup_denoised.nii'])
    unix(['/usr/local/mrtrix3/bin/mrcalc ',dti_tp,'.nii dti7t/dti_topup_denoised.nii -subtract dti7t/dti_topup_residuals.nii'])
     unix('/usr/local/mrtrix3/bin/mrdegibbs -nthreads 14 -axes 0,1 dti7t/dti_topup_denoised.nii dti7t/dti_topup_unring.nii')
%fixB0RingingArtifact('dti7t/dti_topup_denoised.nii','dti_topup_fake','dti7t/dti_topup_unring.nii');

else
     unix('/usr/local/mrtrix3/bin/mrdegibbs -nthreads 14 dti_topup_raw.nii dti7t/dti_topup_unring.nii')
%fixB0RingingArtifact('dti_topup_raw.nii','dti_topup_fake','dti7t/dti_topup_unring.nii');
end; clear cmdout;

%% 
% TOPUP requires even number of slices: check and pad if necessary
% NB: also pads topup file
A = load_untouch_nii('dti7t/dti_unring.nii');
nrslices = A.hdr.dime.dim(4);
if mod(nrslices,2)==1
    A.img(:,:,end+1,:)= zeros(size(A.img,1),size(A.img,2),1,size(A.img,4));
    A.hdr.dime.dim(4)=A.hdr.dime.dim(4)+1;
    save_untouch_nii(A,'dti7t/dti.nii');
    B = load_untouch_nii('dti7t/dti_topup_unring.nii');
    B.img(:,:,end+1,:)= zeros(size(B.img,1),size(B.img,2),1,size(B.img,4));
    B.hdr.dime.dim(4)=B.hdr.dime.dim(4)+1;
    save_untouch_nii(B,'dti7t/dti_tp.nii');
    
    clear A B
else
    copyfile('dti7t/dti_topup_unring.nii','dti7t/dti_tp.nii')
    movefile('dti7t/dti_unring.nii','dti7t/dti.nii')
    clear A
end
delete('dti7t/dti_denoised.nii');

copyfile([dti,'.bval'],'dti7t/dti.bval')
copyfile([dti,'.bvec'],'dti7t/dti.bvec')

clear dti dti_tp nrslices;
cd('dti7t')

%%
% FSL TOPUP and EDDY preprocessing
% NB: before you run it, check the dicoms mentioned in the shell script are
% correct! - they are used to calculate the total readout time


% There are multiple scripts that are very similar to each other. Uncomment
% the one you need and comment the others:

% % % unix(['~/Documents/dti_topup_CPU_only.sh ',dcm,' ',dcm_tp])
%%%unix(['~/Documents/dti_topup_eddy.sh ',dcm,' ',dcm_tp])
!~/Documents/dti_topup_eddy_noROtimecalc.sh


 !mv dti.bvec uncorr_dti.bvec

 !mv eddy_corr_data.eddy_rotated_bvecs dti.bvec

 !mv bet_meanB0_AP_mask.nii.gz nodif_brain.nii.gz
 !mv bet_meanB0_AP_mask_mask.nii.gz nodif_brain_mask.nii.gz
 
 A=load_untouch_nii('eddy_corr_data.nii.gz');
 A.img(A.img(:,:,:,:)<0)=0;
 save_untouch_nii(A,'data.nii.gz')
 clear A
 
%% 
% adaptiveLMMS
%('eddy_corr_data_nonNeg.nii','nodif_brain_mask.nii.gz','dti','filt_'); 

!mv dti.bval bvals
!mv dti.bvec bvecs

% fit diffusion parameters for probabilistic tractography:  FSL's Bedpostx
% and fit a diffusion tensor: FSL's DTIFIT

unix('dtifit -k data.nii -r bvecs -b bvals -m nodif_brain_mask -o dtifit --sse -w')

unix('. /opt/sge/pavlov/common/settings.sh; bedpostx .')