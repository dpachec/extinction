%% imports ACPC coordinates and convert to MNI 

load_paths
mkdir(path_out);

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18','c_sub20'};
       
       
for subji = 1:length(allsubs)
    
    clearvars -except allsubs path_data path_coord path_out subji 
    sel_sub=allsubs{subji};
    info_file=strcat(path_coord,sel_sub,'/', sel_sub, '_datainfo');
    load(info_file)  
    fmri_acpc = ft_read_mri([path_coord sel_sub '/' sel_sub '_T1_fs.nii']);
    fmri_acpc.coordsys = 'acpc';
    elec_acpc_f = datainfo.elec_info.elec_ct_mr; 
    
    % % % Register the subject’s brain to the standard MNI brain ~aprox 1:30min
    cfg = [];
    cfg.nonlinear = 'yes';
    cfg.spmversion = 'spm12';
    cfg.spmmethod  = 'new';
    fmri_mni = ft_volumenormalise(cfg, fmri_acpc);
    

    % % % Use the resulting deformation parameters to obtain the electrode positions in standard MNI space

    elec_mni_frv = elec_acpc_f;
    elec_mni_frv.elecpos = ft_warp_apply(fmri_mni.params, elec_acpc_f.elecpos, 'individual2sn');
    elec_mni_frv.chanpos = ft_warp_apply(fmri_mni.params, elec_acpc_f.chanpos, 'individual2sn');
    elec_mni_frv.coordsys = 'mni';
    
    % % save normalized electrode to file
    save([path_coord sel_sub '/' sel_sub '_elec_mni_frv.mat'], 'elec_mni_frv');
    
end


%% load edf and store only channels in Marie's data

%sub01
cd /Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/raw_data/c_sub01/ieeg
EEG = pop_biosig('DBX-Learning test.edf', 'importevent', 'off'); 
%% plot 
chanids = [49:49];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'winlength', 10, 'spacing', 5000000, 'events', EEG.event);


%%

chans_edf = [{EEG.chanlocs.labels}]';
chans_edf = erase(chans_edf, 'POL '); %deletes over the whole array



[id1 id2 id3] = intersect (chans_edf, elec_mni_frv.label)


%% 
clearvars -except elec_mni_frv
%atlas = ft_read_atlas('/Users/danielpacheco/Documents/MATLAB/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii')
atlas = ft_read_atlas('/Users/danielpacheco/Documents/MATLAB/fieldtrip/template//atlas/afni/TTatlas+tlrc.HEAD')
%atlas = ft_read_atlas('/Users/danielpacheco/Desktop/mni305_lin_nifti/aparc aseg-in-rawavg.nii')
atlas.coordsys = 'mni'; 


for chani = 1:size(elec_mni_frv.elecpos, 1)
    
    cfg = []; 
    cfg.roi = elec_mni_frv.chanpos(chani, :) ; 
    cfg.atlas = atlas; 
    cfg.inputcoord = 'mni'; 
    cfg.output = 'label'; 
    cfg.maxqueryrange = 5; 
    labels = ft_volumelookup(cfg, atlas); 
    find(labels.count)
    [~, indx] = max(labels.count);
    myCoord(chani, :) = labels.name(indx)
    
    
end

















%%