%% Trial level metric of context specificity - SAME EXPERIMENTAL PHASE
%%


clear , clc


listF2sav = {


'trlSTA_AMY_CE_1-44_1_0_500-50';
'trlCTX_AMY_CE_1-44_1_0_500-50';
'trlSTA_TMP_CE_1-44_1_0_500-50';
'trlCTX_TMP_CE_1-44_1_0_500-50';

% 'trlSTA_HPC_CE_1-44_1_0_500-50';
% 'trlCTX_HPC_CE_1-44_1_0_500-50';
% 'trlSTA_PFC_CE_1-44_1_0_500-50';
% 'trlCTX_PFC_CE_1-44_1_0_500-50';
% 'trlSTA_OCC_CE_1-44_1_0_500-50';
% 'trlCTX_OCC_CE_1-44_1_0_500-50';
% 'trlSTA_OFC_CE_1-44_1_0_500-50';
% 'trlCTX_OFC_CE_1-44_1_0_500-50';




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




%% Trial level metric of context specificity - 
% FOR TEST PHASE 

clear , clc


listF2sav = {

'trlSTA_AMY_CET_3-54_1_0_500-100';
'trlSTA_HPC_CET_3-54_1_0_500-100';
'trlSTA_PFC_CET_3-54_1_0_500-100';
'trlSTA_OCC_CET_3-54_1_0_500-100';
'trlSTA_OFC_CET_3-54_1_0_500-100';
'trlSTA_TMP_CET_3-54_1_0_500-100';


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

        case 'trlSTA'
            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear itstaTr
                if ~isempty(nRDM)

                    allTrT = find(ids(:, 2)==3); 
                    firstT = allTrT(1);
                    % loop through every trial during the test
                    count= 1; 
                    for triali = firstT:size(nRDM, 1)
                        currTr = ids(triali, 5); 
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 2 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        itstaTr(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    itstaTRALL{subji, 1} = itstaTr;
                    ids2sav = ids(ids(:, 2) == 3,:); 
                    itstaTRALL{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'itstaTRALL');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% Trial level metric of context specificity 
% FOR EXTINCTION PHASE 

clear , clc


listF2sav = {

'trlSTA_AMY_CAE_3-54_1_0_500-100';
'trlSTA_HPC_CAE_3-54_1_0_500-100';
'trlSTA_PFC_CAE_3-54_1_0_500-100';
'trlSTA_OCC_CAE_3-54_1_0_500-100';
'trlSTA_OFC_CAE_3-54_1_0_500-100';
'trlSTA_TMP_CAE_3-54_1_0_500-100';

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

        case 'trlSTA'
            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear itstaTr
                if ~isempty(nRDM)

                    allTrT = find(ids(:, 2)==2); 
                    firstT = allTrT(1);
                    % loop through every trial during the test
                    count= 1; 
                    for triali = firstT:size(nRDM, 1)
                        currTr = ids(triali, 5); 
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 1 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        itstaTr(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    itstaTRALL{subji, 1} = itstaTr;
                    ids2sav = ids(ids(:, 2) == 2,:); 
                    itstaTRALL{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'itstaTRALL');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))

%% Trial level metric of context specificity 
% FOR ACQUISITION AND TEST 

clear , clc


listF2sav = {

'trlSTA_AMY_CAET_3-54_1_0_500-100';
'trlSTA_HPC_CAET_3-54_1_0_500-100';
'trlSTA_PFC_CAET_3-54_1_0_500-100';
'trlSTA_OCC_CAET_3-54_1_0_500-100';
'trlSTA_OFC_CAET_3-54_1_0_500-100';
'trlSTA_TMP_CAET_3-54_1_0_500-100';

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

        case 'trlSTA'
            for subji = 1:length(nRDMs)
                nRDM = nRDMs{subji,1}; 
                ids = double(string(nRDMs{subji,2})); 
                clear itstaTr
                if ~isempty(nRDM)

                    allTrT = find(ids(:, 2)==3); 
                    firstT = allTrT(1);
                    % loop through every trial during the test
                    count= 1; 
                    for triali = firstT:size(nRDM, 1)
                        currTr = ids(triali, 5); 
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 1 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        itstaTr(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    itstaTRALL{subji, 1} = itstaTr;
                    ids2sav = ids(ids(:, 2) == 3,:); 
                    itstaTRALL{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'itstaTRALL');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))
%%

