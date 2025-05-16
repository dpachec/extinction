%% Export power
% Load traces (path specified in load_paths_EXT.m) and export power
% Export power locked to video onset (_V_) or cue onset (_C_)
clear , clc

listF2sav = {


'POW_AMY_C_10';
'POW_TMP_C_10';
'POW_PFC_C_10';
'POW_OFC_C_10';
'POW_AMY_C_10';


'POW_HPC_V_10';
'POW_TMP_V_10';
'POW_PFC_V_10';
'POW_OFC_V_10';
'POW_AMY_V_10';

};   

paths = load_paths_EXT; 
t1 = datetime; 

for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
            
    f2sav       = listF2sav{listi}; 
    cfg         = getParams_EXT(f2sav);

    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '    ' listF2sav{listi} '   Subject > ' num2str(subji)]);
        EEG = ALLEEG{subji};
        if ~isempty(EEG)
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); 

            if strcmp(cfg.roi, 'aHPC')
                EEG = select_ant_hpc_chans(EEG); 
            end

            if ~isempty(EEG.data)
                Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
                Ev2 = cat(1, Ev1{:});
                
                EEG = extract_power_EXT(EEG, cfg); 
                EEG = normalize_EXT(EEG);  %across trials
                oneListPow = EEG.power(:, :, : ,251:550); 
                POW{subji, 1} = oneListPow;
                POW{subji, 2} = Ev2; 
            end
            
        end
        
    end
   
    save([ paths.results.POWfromRT f2sav '.mat'], 'POW');
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end












%%


























%%


