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



%% produce meshes from AVERAGE MNI .nii in freesurver


clear

fshome = '/Applications/freesurfer/7.3.2';
subdir = ['/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/brainMNI/'];
mrfile = ['/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/brainMNI/anat/mni_icbm152_t1_tal_nlin_sym_09b_hires.nii'];
%this brain is downloaded from https://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009
% take the non-linear symmetric 2009 b model, as in the allan instititu
% converting this in freesurfer gives the same surfaces tony had in bx3
% (previous nat comm and pnas paper use this

system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/anat/tmp.nii'] '; ' ...2
'recon-all -i ' [subdir '/anat/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all']);






%% convert to OBJ
clear 

fshome = '/Applications/freesurfer/7.3.2';
pialFile = '/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/brainMNI/freesurfer/surf/rh.pial';
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mris_convert ' pialFile ' /Users/danielpacheco/Documents/iEEG_data_analysis/extinction/brainMNI/freesurfer/surf/rh.pial.stl; ']);

%taking rh.pial or rh.pial.T1 gives identical results (same model)
% stl files can be read in matlab and easily converted to obj fbx etc
% importing stl in maya requires loading the plug in windows,
% settings/preferences, plug-in manager, stlTranslator.mll







%% plot 

atlas = ft_read_atlas('/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/brainMNI/freesurfer/mri/aparc+aseg.mgz');
atlas.coordsys = 'mni';
cfg            = [];
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
cfg.roi        = {'Left-Hippocampus'};
mask_rha     = ft_volumelookup(cfg, atlas);
seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha_subject = ft_prepare_mesh(cfg, seg);


%%
figure(2)
ft_plot_mesh(mesh_rha_subject,  'facealpha', .5);
%ft_plot_sens(elec_acpc_f);
title('Hippocampus MNI coordinates');



















