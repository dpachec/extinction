%% 
clear, close all
paths = load_paths; 
currentPath = pwd; 
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub29' };

%% epoch traces


for subji = 1:length(allsubs)
    
    sub = allsubs{subji}; 
    cd([ paths.ds '/' sub])
    log_list = dir('*mat'); log_list = {log_list.name};
    load ([sub '_downSampiEEG.mat']);
    
    Ev = [EEG.event.type]'; 
    Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
    Ev2 = cat(1, Ev1{:})
    ids = strcmp(Ev2(:, 10), 'T')
    EEG.event = EEG.event(ids)
    EEG_b = pop_epoch( EEG, {}, [-4 8], 'newname', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
    EEG = rem_EEGLAB_fields(EEG_b);
    





end



%% 
d2p = squeeze(EEG.data(1, :, 1)); 
figure;plot(d2p)











%% 