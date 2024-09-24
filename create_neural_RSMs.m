%%
%% CREATE NEURAL RSMs
%nRSM_Region_CUE-TaskPeriod_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_
clear , clc

listF2sav = {


'nRSM_aHPC_CAT_1-44_1_0_500-50';
% 'nRSM_HPC_CA_1-44_1_0_500-50';
% 'nRSM_OFC_CA_1-44_1_0_500-50';
% 'nRSM_OCC_CA_1-44_1_0_500-50';
% 'nRSM_TMP_CA_1-44_1_0_500-50';
% 'nRSM_PFC_CA_1-44_1_0_500-50';
% 'nRSM_AMY_CA_1-44_1_0_500-50';




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
            if strcmp(cfg.period, 'A')
                ids2rem = ids(:, 2)~=1; 
            elseif strcmp(cfg.period, 'E')
                ids2rem = ids(:, 2)~=2; 
            elseif strcmp(cfg.period, 'T')
                ids2rem = ids(:, 2)~=3; 
            elseif strcmp(cfg.period, 'ET')
                ids2rem = ids(:, 2)==1; 
            elseif strcmp(cfg.period, 'AE')
                ids2rem = ids(:, 2)==3; 
            elseif strcmp(cfg.period, 'AT')
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
