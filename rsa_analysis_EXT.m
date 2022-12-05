%% RSA ANALYSIS
%% Load data

clear 
paths = load_paths_EXT; 
file2load = ['allS_' 'Amygdala' '_C'];

load ([paths.results.power file2load]); 


%% compute model fit
%freqs_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_win-width_mf
clearvars -except ALLEEG paths file2load

f2sav = '30-54_0_0_0_20-1'; 
cfg = getParams_EXT(f2sav);


for subji = 1:length(ALLEEG)
    
    subji

    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        if ndims(EEG.power) == 3
            tmph(:, 1, :, :)  = EEG.power; 
            EEG.power = tmph; 
        end

        ids1 =  strcmp(Ev2(:, 2), '1') & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') ;
        %ids1 =  strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1')  ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3' )  | strcmp(Ev2(:, 6), '2') ) ;
        allCSPIts = EEG.power(ids1, :, : ,201:500); 
        allCSMIts = EEG.power(ids2, :, : ,201:500); 
        
        [allCSPIts] = remove_nans_EXT(allCSPIts);
        [allCSMIts] = remove_nans_EXT(allCSMIts);

        allIts = cat(1, allCSPIts, allCSMIts);
        
        neuralRDM = createNeuralRDMs_EXT(cfg, allIts);
        modelRDM  = createModelRDM_EXT(neuralRDM, allCSPIts);
        
        nnFit{subji,1}              = fitModel_WM(neuralRDM, modelRDM, cfg.fitMode); 
        %nnFit{subji,2}              = cfg_contrasts.oneListIds_c; 

        
    end

end


filename = [paths.results.rsa 'CM_' file2load];
save(filename, 'nnFit');


cd (paths.github)

disp('done');

%% plot grand average fit

d2p = cell2mat(nnFit);
md2p = mean(d2p);
std2p = std(d2p, [], 1); 
set2p = std2p / sqrt(size(d2p, 1));
[h p ci ts] = ttest(d2p)
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 



%plot(md2p); hold on; 
%plot(h)



times = 1:size(d2p,2); 
shadedErrorBar(times, md2p, set2p, 'r', 1); hold on; 
plot (times, hb, 'Linewidth', 4)
set(gca, 'xtick', [1 90 240], 'xticklabels', {'-1' '0' '1.5'}, 'xlim', [51 251])
set(gca, 'FontSize', 12);
plot([90 90],get(gca,'ylim'), 'k','lineWidth',1, 'Color', [.5 .5 .5]);
plot(get(gca,'xlim'), [0 0 ], 'k','lineWidth',1, 'Color', [.5 .5 .5]);










%%




























%%