%%
%% CREATE NEURAL RSMs
%nRSM_Region_CUE-TaskPeriod_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_
clear , clc

listF2sav = {

'nRSM_aHPC_CT_3-54_1_0_500-100';
'nRSM_HPC_CT_3-54_1_0_500-100';
'nRSM_OFC_CT_3-54_1_0_500-100';
'nRSM_OCC_CT_3-54_1_0_500-100';
'nRSM_TMP_CT_3-54_1_0_500-100';
'nRSM_PFC_CT_3-54_1_0_500-100';
'nRSM_AMY_CT_3-54_1_0_500-100';


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
            if cfg.period == 'A'
                ids2rem = ids(:, 2)~=1; 
            elseif cfg.period == 'E'
                ids2rem = ids(:, 2)~=2; 
            elseif cfg.period == 'T'
                ids2rem = ids(:, 2)~=3; 
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
