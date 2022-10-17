%%
clear
paths = load_paths; 

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18','c_sub20'};


for subji = 3:length(allsubs)

subjID = char(allsubs(subji));
%mri = ft_read_mri([path_coord subjID '/' subjID '_ct_acpc_f.nii']);
fshome = '/Applications/freesurfer/7.3.2';
subdir = [paths.coord subjID '/'];
%mrfile = [path_coord subjID '/' subjID '_T1_fs.nii'];
mrfile = [paths.coord subjID '/anat/' subjID '_MR_acpc.nii'];
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/anat/tmp.nii'] '; ' ...2
'recon-all -i ' [subdir '/anat/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all']);

end


%%
