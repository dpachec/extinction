%% Trial level metric of context specificity - SAME EXPERIMENTAL PHASE
%%


clear , clc


listF2sav = {

'trlSTA_TMP_VE_1-44_1_0_500-50';
'trlCTX_TMP_VE_1-44_1_0_500-50';

% 'trlSTA_HPC_VE_1-44_1_0_500-50';
% 'trlCTX_HPC_VE_1-44_1_0_500-50';
% 'trlSTA_PFC_VE_1-44_1_0_500-50';
% 'trlCTX_PFC_VE_1-44_1_0_500-50';
% 'trlSTA_OCC_VE_1-44_1_0_500-50';
% 'trlCTX_OCC_VE_1-44_1_0_500-50';
% 'trlSTA_OFC_VE_1-44_1_0_500-50';
% 'trlCTX_OFC_VE_1-44_1_0_500-50';
% 'trlSTA_AMY_VE_1-44_1_0_500-50';
% 'trlCTX_AMY_VE_1-44_1_0_500-50';


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




%% Trial level metric of item stability between extinction and test > loop through every trial during extinction

clear , clc


listF2sav = {

'trlSTA_aHPC_CET_1-44_1_0_500-50';

% 'trlSTA_AMY_CET_1-44_1_0_500-50';
% 'trlSTA_HPC_CET_1-44_1_0_500-50';
% 'trlSTA_PFC_CET_1-44_1_0_500-50';
% 'trlSTA_OCC_CET_1-44_1_0_500-50';
% 'trlSTA_OFC_CET_1-44_1_0_500-50';
% 'trlSTA_TMP_CET_1-44_1_0_500-50';
% 

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
                clear reinstET
                if ~isempty(nRDM)

                    allTrE = find(ids(:, 2)==2); 
                    
                    % loop through every trial during extinction
                    count= 1; 
                    for triali = 1:length(allTrE)
                        t2u = allTrE(triali);
                        currTr = ids(t2u, 5); 
                        %identify the same trial during test
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 3 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        reinstET(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    trlSTA_ET{subji, 1} = reinstET;
                    ids2sav = ids(ids(:, 2) == 2,:); 
                    trlSTA_ET{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'trlSTA_ET');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))

%% Trial level metric of item stability between extinction and test > loop through every trial during test

clear , clc


listF2sav = { %note that in this case the names are modified at the end CET > CTE

'trlSTA_AMY_CET_1-44_1_0_500-50';
'trlSTA_HPC_CET_1-44_1_0_500-50';
'trlSTA_PFC_CET_1-44_1_0_500-50';
'trlSTA_OCC_CET_1-44_1_0_500-50';
'trlSTA_OFC_CET_1-44_1_0_500-50';
'trlSTA_TMP_CET_1-44_1_0_500-50';


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
                clear reinstTE
                if ~isempty(nRDM)

                    allTrT = find(ids(:, 2)==3); 
                    
                    % loop through every trial during extinction
                    count= 1; 
                    for triali = 1:length(allTrT)
                        t2u = allTrT(triali);
                        currTr = ids(t2u, 5); 
                        %identify the same trial during extinction
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 2 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        reinstTE(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    trlSTA_TE{subji, 1} = reinstTE;
                    ids2sav = ids(ids(:, 2) == 3,:); 
                    trlSTA_TE{subji, 2} = ids2sav;
                end
            end

            % change name and save as CTE
            f2sav = strrep(f2sav, 'CET', 'CTE'); 
            save([ paths.results.trial_based f2sav '.mat'], 'trlSTA_TE');
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))



%% Trial level metric of item stability between acquisition and test > loop through every trial during acquisition

clear , clc


listF2sav = {

'trlSTA_AMY_CAT_1-44_1_0_500-50';
'trlSTA_HPC_CAT_1-44_1_0_500-50';
'trlSTA_aHPC_CAT_1-44_1_0_500-50';
'trlSTA_PFC_CAT_1-44_1_0_500-50';
'trlSTA_OCC_CAT_1-44_1_0_500-50';
'trlSTA_OFC_CAT_1-44_1_0_500-50';
'trlSTA_TMP_CAT_1-44_1_0_500-50';


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
                clear reinstET
                if ~isempty(nRDM)

                    allTrA = find(ids(:, 2)==1); 
                    
                    % loop through every trial during extinction
                    count= 1; 
                    for triali = 1:length(allTrA)
                        t2u = allTrA(triali);
                        currTr = ids(t2u, 5); 
                        %identify the same trial during test
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 3 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        reinstAT(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    trlSTA_AT{subji, 1} = reinstAT;
                    ids2sav = ids(ids(:, 2) == 1,:); 
                    trlSTA_AT{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'trlSTA_AT');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% Trial level metric of item stability between acquisition and test > IMPORTANT !! %note that in this case the names are modified at the end CAT > CTA

clear , clc


listF2sav = { %note that in this case the names are modified at the end CAT > CTA

'trlSTA_AMY_CAT_1-44_1_0_500-50';
'trlSTA_HPC_CAT_1-44_1_0_500-50';
'trlSTA_aHPC_CAT_1-44_1_0_500-50';
'trlSTA_PFC_CAT_1-44_1_0_500-50';
'trlSTA_OCC_CAT_1-44_1_0_500-50';
'trlSTA_OFC_CAT_1-44_1_0_500-50';
'trlSTA_TMP_CAT_1-44_1_0_500-50';


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
                clear reinstAT
                if ~isempty(nRDM)

                    allTrT = find(ids(:, 2)==3); 
                    
                    % loop through every trial during test
                    count= 1; 
                    for triali = 1:length(allTrT)
                        t2u = allTrT(triali);
                        currTr = ids(t2u, 5); 
                        %identify the same trial during acquisition
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 1 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        reinstAT(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    trlSTA_TA{subji, 1} = reinstAT;
                    ids2sav = ids(ids(:, 2) == 3,:); 
                    trlSTA_TA{subji, 2} = ids2sav;
                end
            end

            f2sav = strrep(f2sav, 'CAT', 'CTA'); 
            save([ paths.results.trial_based f2sav '.mat'], 'trlSTA_TA');
    
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

'trlSTA_AMY_CAT_1-44_1_0_500-50';
'trlSTA_HPC_CAT_1-44_1_0_500-50';
'trlSTA_PFC_CAT_1-44_1_0_500-50';
'trlSTA_OCC_CAT_1-44_1_0_500-50';
'trlSTA_OFC_CAT_1-44_1_0_500-50';
'trlSTA_TMP_CAT_1-44_1_0_500-50';

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

                    allTrA = find(ids(:, 2)==1); 
                    % loop through every trial during the test
                    count= 1; 
                   for triali = 1:length(allTrA)
                        t2u = allTrA(triali);
                        currTr = ids(t2u, 5); 
                        sameITIDs = (ids(:, 5) == currTr) & ( ids(:, 2) == 3 ) ;
                        % sameITIDs(triali) = 0; % no need to exclude the same trial here (diff exp phases)
                        itstaTr(count, :,:) = squeeze( mean(nRDM(triali, sameITIDs,:,:), 2) );
                        count = count+1; 
                    end
                    trlSTA_AT{subji, 1} = itstaTr;
                    ids2sav = ids(ids(:, 2) == 3,:); 
                    trlSTA_AT{subji, 2} = ids2sav;
                end
            end
            save([ paths.results.trial_based f2sav '.mat'], 'trlSTA_AT');            
    
    end
            
end



t2 = datetime; 
etime(datevec(t2), datevec(t1))
%%

