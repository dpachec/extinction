%% Trial level metric of context specificity
%%


clear , clc


listF2sav = {


'trlSTA_AMY_CE_13-29_1_0_500-100';
'trlCTX_AMY_CE_13-29_1_0_500-100';
'trlSTA_HPC_CE_13-29_1_0_500-100';
'trlCTX_HPC_CE_13-29_1_0_500-100';
'trlSTA_PFC_CE_13-29_1_0_500-100';
'trlCTX_PFC_CE_13-29_1_0_500-100';
'trlSTA_OCC_CE_13-29_1_0_500-100';
'trlCTX_OCC_CE_13-29_1_0_500-100';
'trlSTA_OFC_CE_13-29_1_0_500-100';
'trlCTX_OFC_CE_13-29_1_0_500-100';
'trlSTA_TMP_CE_13-29_1_0_500-100';
'trlCTX_TMP_CE_13-29_1_0_500-100';




};   


t1 = datetime; 
for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    f2t = strsplit(f2sav, '_'); 
    paths = load_paths_EXT; 
    
    % load neural RSM
    load ([paths.results.nRDMs 'nRSM' listF2sav{listi}(7:end)])


    switch f2t{1}
        case 'trlCTX'

            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear ctxTr
                if ~isempty(nRDM)
                    for triali = 1:size(nRDM, 1)
                        currTr = ids(triali, 3); 
                        sameCTXIDs = ids(:, 3) == currTr;
                        sameCTXIDs(triali) = 0; %exclude the same trial by definition 1
                        diffCTXIDs = ids(:, 3) ~= currTr;
                        ctxTr(triali, :,:) = squeeze( mean(nRDM(triali, sameCTXIDs,:,:), 2) - mean(nRDM(triali, diffCTXIDs,:,:), 2) );
                    end
                    ctxTRALL{subji, 1} = ctxTr;
                    ctxTRALL{subji, 2} = ids; 
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'ctxTRALL');

        case 'trlSPE'
            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear itspecTr
                if ~isempty(nRDM)
                    for triali = 1:size(nRDM, 1)
                        currTr = ids(triali, 5); 
                        sameCTXIDs = ids(:, 5) == currTr;
                        sameCTXIDs(triali) = 0; %exclude the same trial by definition 1
                        diffCTXIDs = ids(:, 5) ~= currTr;
                        itspecTr(triali, :,:) = squeeze( mean(nRDM(triali, sameCTXIDs,:,:), 2) - mean(nRDM(triali, diffCTXIDs,:,:), 2) );
                    end
                    itspecTRALL{subji, 1} = itspecTr;
                    itspecTRALL{subji, 2} = ids; 
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'itspecTRALL');

        case 'trlSTA'
            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear itstaTr
                if ~isempty(nRDM)
                    for triali = 1:size(nRDM, 1)
                        currTr = ids(triali, 5); 
                        sameCTXIDs = ids(:, 5) == currTr;
                        sameCTXIDs(triali) = 0; %exclude the same trial by definition 1
                        itstaTr(triali, :,:) = squeeze( mean(nRDM(triali, sameCTXIDs,:,:), 2) );
                    end
                    itstaTRALL{subji, 1} = itstaTr;
                    itstaTRALL{subji, 2} = ids;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'itstaTRALL');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))




%%

