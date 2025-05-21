%%
%% CREATE NEURAL RSMs at different task stages
% Code to generate the file 'nRSM'_ROI_Locked2_freqs2include_avTimeFeatVect_freqResolv(0-1)_win-width_mf_"
% nRSM : identifier of the file
% ROI = HPC, TMP, OFC, PFC, AMY
% Locked 2: CUE-TASK PERIOD. E.g.; VE = Video during extinction; CE: Cue during extinction; VA: Video during acquisition; CA: cue during acquisition
% freqs2include = 1-44 (all frequencies in the 1-150Hz; see extract_power_EXT.m)
% avTimeFeatVect: 0: all time points are included in the representational feature vectors); 1: average power values within each window)
% freqREsolved: 0: all frequencies are included in the representational feature vectors); 1: feature vectors are extracted for each frequency
% win-width: window length (datapoints)
% mf: sliding in (datapoints)  
% see get_params.m
% this script loads power data extracted in with the create_POW.m script. These files can also be downloaded from https://osf.io/xfprt/POWfromRT/
clear , clc

listF2sav = {

 'nRSM_HPC_VE_1-44_1_0_500-50';
 'nRSM_OFC_VE_1-44_1_0_500-50';
 'nRSM_TMP_VE_1-44_1_0_500-50';
 'nRSM_PFC_VE_1-44_1_0_500-50';
 'nRSM_AMY_VE_1-44_1_0_500-50';




};   

paths = load_paths_EXT; 
t1 = datetime; 

for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
        
    f2sav       = listF2sav{listi}; 
    cfg         = getParams_EXT(f2sav);
    
    
    ALLPOW = load ([paths.results.POWfromRT cfg.powF2load]); 
     
    
    for subji = 1:length(ALLPOW.POW)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        POW = ALLPOW.POW{subji, 1};
        ids = double(string(ALLPOW.POW{subji, 2}));
        

        if ~isempty(POW)
            % % select experimental phase
            if strcmp(cfg.period, 'A') % acquisition
                ids2rem = ids(:, 2)~=1; 
            elseif strcmp(cfg.period, 'E') % extinction 
                ids2rem = ids(:, 2)~=2; 
            elseif strcmp(cfg.period, 'T') % test
                ids2rem = ids(:, 2)~=3; 
            elseif strcmp(cfg.period, 'ET') % extinction to test
                ids2rem = ids(:, 2)==1; 
            elseif strcmp(cfg.period, 'AE') % acquisition to extinction
                ids2rem = ids(:, 2)==3; 
            elseif strcmp(cfg.period, 'AT') % acquisitoin to test
                ids2rem = ids(:, 2)==2; 
            elseif strcmp(cfg.period, 'AET')
                ids2rem = []; 
            end


            ids(ids2rem,:) = [];
            POW(ids2rem,:,:,:) = []; 
            cfg.POW = POW; 
            nRDMs{subji, 1} = createNeuralRDMs_EXT(cfg);
            nRDMs{subji, 2} = ids; 
        end
    end

    save([ paths.results.nRDMs f2sav '.mat'], 'nRDMs');
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end



%%
