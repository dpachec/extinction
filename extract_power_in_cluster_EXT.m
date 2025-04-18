%% % % % % % % % Extract amygdala THETA power in the cluster
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 

load ([paths.results.power file2load]); 


%% EXTRACT SINGLE TRIAL POWER LEVELS

clearvars -except ALLEEG paths clustinfo nL 
close all, clc

paths = load_paths_EXT; 

%load ([paths.results.clusters 'AMY_THETA_CLUST_px4'])
%load ([paths.results.clusters 'AMY_POWER_CLUST_px42'])
load _44_clustinfo_AMY_THETA_px32.mat
px = 32; 

% % % % %  check the cluster 
% d2p = zeros(44, 200); 
% d2p(clustinfo.PixelIdxList{px}) = 1; 
% contour(d2p)

for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    clear thPow
    if ~isempty(EEG) 
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2');
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        %length(find(ids))
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            
            % % % Theta only 
            cTR = squeeze(mean(powH(triali, :, 1:44,276:475), 2)); %mean across channels
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{px}), 'all');



            % % % All frequencies
            %cTR = squeeze(mean(powH(triali, :, 3:54, 251:500), 2));
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{42}), 'all');


            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            %cTR = squeeze(mean(powH(triali, :, 3:8, 441:460), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
            nNan{subji,:} = sum(isnan(thPow));

           allPOWAMY{subji, 1} = thPow; 
           ids2sav = Ev2(ids, :); 
           allPOWAMY{subji, 2} = ids2sav; 


    end

    

end


save ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'], 'allPOWAMY', 'nNan'); 






%%
















