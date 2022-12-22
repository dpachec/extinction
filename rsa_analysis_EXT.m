%% RSA ANALYSIS
%% Load data

clear 
paths = load_paths_EXT; 
file2load = ['allS_' 'Hippocampus' '_C'];
%file2load = 'allS_inferiortemporal_middletemporal_superiortemporal_bankssts_ctx-lh-fusiform_ctx-lh-temporalpole_inferiorparietal_lateraloccipital_lingual_parahippocampal_cuneus_pericalcarine_C.mat'
load ([paths.results.power file2load]); 


%% compute neural RDM 
%freqs_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_win-width_mf
clearvars -except ALLEEG paths file2load

f2sav = '39-54_1_0_0_50-1'; 
cfg = getParams_EXT(f2sav);


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        if ndims(EEG.power) == 3
            tmph(:, 1, :, :)  = EEG.power; 
            EEG.power = tmph; 
        else %pick just the deepest amygdala electrode
            
            %nChans = size(EEG.power, 2); 
            %if nChans > 3
            %    EEG.power = EEG.power(:, 1:3, :, :);
            %end


        end
        

        ids1 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '1');
        ids2 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '2');
        ids3 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        all1 = EEG.power(ids1, :, : ,201:550); 
        all2 = EEG.power(ids2, :, : ,201:550); 
        all3 = EEG.power(ids3, :, : ,201:550); 
        
        [all1] = remove_nans_EXT(all1);
        [all2] = remove_nans_EXT(all2);
        [all3] = remove_nans_EXT(all3);
        
        allIts = cat(1, all1, all2, all3);

        minTR = 3;  
        if size(all3, 1) > minTR  % CS- condition has fewer trials
        
            %neuralRDM1 = createNeuralRDMs_EXT(cfg, all1);
            %neuralRDM2 = createNeuralRDMs_EXT(cfg, all2);
            %neuralRDM3 = createNeuralRDMs_EXT(cfg, all3);
            
            neuralRDMALL = createNeuralRDMs_EXT(cfg, allIts);
            
            %neuralRDM1 = remDiagRDM_TS(neuralRDM1); 
            %neuralRDM2 = remDiagRDM_TS(neuralRDM2); 
            %neuralRDM3 = remDiagRDM_TS(neuralRDM3); 
            neuralRDMALL = remDiagRDM_TS(neuralRDMALL); 
    
    
            sRDM = size(neuralRDMALL, 1); 
            sRDM1= size(all1, 1);
            sRDM2= size(all2, 1);
            
            % % % CHECK THAT TEMPLATE USED ARE CORRECT 
% % %             SI1_T = zeros(sRDM);SI2_T = zeros(sRDM);SI3_T = zeros(sRDM);
% % %             SI1_T(1:sRDM1,1:sRDM1) = 1; 
% % %             SI2_T(sRDM1+1:sRDM2+sRDM1,sRDM1+1:sRDM2+sRDM1) = 1; 
% % %             SI3_T(sRDM1+sRDM2+1:end,sRDM1+sRDM2+1:end) = 1; 
% % %             
% % %             DISV = zeros(sRDM); 
% % %             DISV(1:sRDM1, sRDM1+1:sRDM2+sRDM1) = 1; 
% % %             DISV(eye(size(DIDV, 1))== 1) = 2
% % %             figure; imagesc(DISV); axis square
% % %             DIDV = zeros(sRDM); 
% % %             DIDV([1:sRDM1 sRDM1+1:sRDM2+sRDM1], [sRDM1+sRDM2+1:end sRDM1+sRDM2+1:end]) = 1; 
% % %             DIDV(eye(size(DIDV, 1))== 1) = 2
% % %             figure; imagesc(DIDV); axis square


            SI1 = neuralRDMALL(1:sRDM1,1:sRDM1, :); 
            SI2 = neuralRDMALL(sRDM1+1:sRDM2+sRDM1,sRDM1+1:sRDM2+sRDM1, :); 
            SI3 = neuralRDMALL(sRDM1+sRDM2+1:end,sRDM1+sRDM2+1:end, :); 
            DISV = neuralRDMALL (1:sRDM1, sRDM1+1:sRDM2+sRDM1,:); 
            DIDV = neuralRDMALL ([1:sRDM1 sRDM1+1:sRDM2+sRDM1], [sRDM1+sRDM2+1:end sRDM1+sRDM2+1:end],:); 


            avCorrSI1(subji, :) = mean(mean(SI1, 1, 'omitnan'), 2, 'omitnan');
            avCorrSI2(subji, :) = mean(mean(SI2, 1, 'omitnan'), 2, 'omitnan');
            avCorrSI3(subji, :) = mean(mean(SI3, 1, 'omitnan'), 2, 'omitnan');
            avCorrDISV(subji, :) = mean(mean(DISV, 1, 'omitnan'), 2, 'omitnan');
            avCorrDIDV(subji, :) = mean(mean(DIDV, 1, 'omitnan'), 2, 'omitnan');


        end
        
    end

end


avCorrSI1 =  avCorrSI1 (any(avCorrSI1 ,2),:);
avCorrSI2 =  avCorrSI2 (any(avCorrSI2 ,2),:);
avCorrSI3 =  avCorrSI3 (any(avCorrSI3 ,2),:);
avCorrDISV =  avCorrDISV (any(avCorrDISV ,2),:);
avCorrDIDV =  avCorrDIDV (any(avCorrDIDV ,2),:);


cd (paths.github)

disp('done');

%% plot 2 conditions

sub2exc = []

avSI1 = avCorrSI1; avSI2 = avCorrSI2; 
avSICP = squeeze(mean(cat(3, avCorrSI1, avCorrSI2), 3)); 
avSICM = avCorrSI3; 
avDISV = avCorrDISV; avDIDV = avCorrDIDV;
avSI1(sub2exc, :) = []; avSI2(sub2exc, :) = []; avSICM(sub2exc, :) = []; avSICP(sub2exc, :) = []; 
avDISV(sub2exc, :) = []; avDIDV(sub2exc, :) = []; 

mavSICP = mean(avSICP); 
mavSICM = mean(avSICM); 
mavDISV = mean(avDISV); 
mavDIDV = mean(avDIDV); 


%times = (-1:.01:2) + .25
times = 1:size(mavSICP,2); 
figure(); 
plot(times, mavSICP, 'r', LineWidth=2); hold on; 
plot(times, mavSICM, 'k', LineWidth=2); hold on; 
%set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)
legend({'CS+' 'CS-'})

figure(); 
plot(times, mavDISV, 'r', LineWidth=2); hold on; 
plot(times, mavDIDV, 'k', LineWidth=2); hold on; 
%set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)
legend({'DISV' 'DIDV' 'CS-'})


%% 

mSICP = mean(avSICP);
stdSICP = std(avSICP, [], 1); 
seSICP = stdSICP / sqrt(size(avSICP, 1));

mSICM = mean(avSICM);
stdSICM = std(avSICM, [], 1); 
seSICM = stdSICM / sqrt(size(avSICM, 1));

diffC = avSICP - avSICM; 
[h p ci ts] = ttest(diffC); 
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 


%times = 1:size(avC1,2); 
times = (-1:.01:2) + .25
shadedErrorBar(times, mSICP, seSICP, 'r', 1); hold on; 
shadedErrorBar(times, mSICM, seSICM, 'k', 1); hold on; 
plot (times, hb, 'Linewidth', 4)
set(gca, 'FontSize', 18);
%set(gca, 'xtick', [1 90 240], 'xticklabels', {'-1' '0' '1.5'}, 'xlim', [51 251])
%plot([90 90],get(gca,'ylim'), 'k','lineWidth',1, 'Color', [.5 .5 .5]);
%plot(get(gca,'xlim'), [0 0 ], 'k','lineWidth',1, 'Color', [.5 .5 .5]);

%% plot 3 items separately

sub2exc = [1]

avC1 = avCorrSI1; avC2 = avCorrSI2; avC3 = avCorrSI3; 
avC1(sub2exc, :) = []; avC2(sub2exc, :) = []; avC3(sub2exc, :) = []; 

mav1 = mean(avC1); 
mav2 = mean(avC2); 
mav3 = mean(avC3); 



% figure(); 
% plot(mCPM - mCPP, 'r', LineWidth=2); hold on; 
% plot(mCPM - mCMM, 'b', LineWidth=2); hold on; 

times = (-1:.01:2) + .25
%times = 1:size(d2p1,2); 
figure(); 
plot(times, mav1, 'r', LineWidth=2); hold on; 
plot(times, mav2, 'r', LineWidth=2); hold on; 
plot(times, mav3, 'k', LineWidth=2); hold on; 
set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)

legend({'CS+1' 'CS+2' 'CS-'})


%% 
md2p1 = mean(avC1);
std2p1 = std(avC1, [], 1); 
set2p1 = std2p1 / sqrt(size(avC1, 1));
md2p2 = mean(avC2);
std2p2= std(avC2, [], 1); 
set2p2 = std2p2 / sqrt(size(avC2, 1));
md2p3 = mean(avC3);
std2p3 = std(avC3, [], 1); 
set2p3 = std2p3 / sqrt(size(avC3, 1));



%times = 1:size(d2p1,2); 
times = (-1:.01:2) + .25
shadedErrorBar(times, md2p1, set2p1, 'r', 1); hold on; 
shadedErrorBar(times, md2p2, set2p2, 'r', 1); hold on; 
shadedErrorBar(times, md2p3, set2p3, 'k', 1); hold on; 
%plot (times, hb, 'Linewidth', 4)
%set(gca, 'xtick', [1 90 240], 'xticklabels', {'-1' '0' '1.5'}, 'xlim', [51 251])
%set(gca, 'FontSize', 12);
%plot([90 90],get(gca,'ylim'), 'k','lineWidth',1, 'Color', [.5 .5 .5]);
%plot(get(gca,'xlim'), [0 0 ], 'k','lineWidth',1, 'Color', [.5 .5 .5]);




%% Analysis 2: across items
%freqs_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_win-width_mf
clearvars -except ALLEEG paths file2load

f2sav = '3-54_0_0_0_50-1'; 
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
        else %pick just the deepest amygdala electrode
            
            %nChans = size(EEG.power, 2); 
            %if nChans > 3
            %    EEG.power = EEG.power(:, 1:3, :, :);
            %end


        end
        

        ids1 =  strcmp(Ev2(:, 2), '1') & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) <= 40;

        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') ;
        %ids2 = strcmp(Ev2(:, 6), '3')  & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) <= 40;
        
        %ids1 =  strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1')  ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3' )  | strcmp(Ev2(:, 6), '2') ) ;
        allCSPIts = EEG.power(ids1, :, : ,201:550); 
        allCSMIts = EEG.power(ids2, :, : ,201:550); 
        
        [allCSPIts] = remove_nans_EXT(allCSPIts);
        [allCSMIts] = remove_nans_EXT(allCSMIts);

        allIts = cat(1, allCSPIts, allCSMIts);
        
        if ~isempty(allIts)

            neuralRDM = createNeuralRDMs_EXT(cfg, allIts);
            
            
            sRDM = size(neuralRDM, 1); 
            nPPits = size(allCSPIts, 1);
            cPM = neuralRDM (1:nPPits,nPPits+1:sRDM,:); 
            cPP = neuralRDM (1:nPPits,1:nPPits,:); 
            cMM = neuralRDM (nPPits+1:sRDM,nPPits+1:sRDM,:); 
    
            % remove diagonal at each time point
            for timei = 1:size(cPP, 3)
                cPPt = cPP(:, :, timei); 
                cPPt(eye(size(cPPt))==1) = nan;
                cPP(:, :, timei) = cPPt; 
                cMMt = cMM(:, :, timei); 
                cMMt(eye(size(cMMt))==1) = nan;
                cMM(:, :, timei) = cMMt; 
                
            end
            avCorrPM(subji, :) = mean(mean(cPM, 1, 'omitnan'), 2, 'omitnan');
            avCorrPP(subji, :) = mean(mean(cPP, 1, 'omitnan'), 2, 'omitnan');
            avCorrMM(subji, :) = mean(mean(cMM, 1, 'omitnan'), 2, 'omitnan');

        end


        
    end

end


avCorrPM =  avCorrPM (any(avCorrPM ,2),:);
avCorrPP =  avCorrPP (any(avCorrPP ,2),:);
avCorrMM =  avCorrMM (any(avCorrMM ,2),:);


cd (paths.github)

disp('done');



%%

mCPM = mean(avCorrPM); 
mCPP = mean(avCorrPP); 
mCMM = mean(avCorrMM); 

% figure(); 
% plot(mCPM - mCPP, 'r', LineWidth=2); hold on; 
% plot(mCPM - mCMM, 'b', LineWidth=2); hold on; 

times = (-1:.01:2) + .25
figure(); 
plot(times, mCPP, 'r', LineWidth=2); hold on; 
plot(times, mCMM, 'b', LineWidth=2); hold on; 
plot(times, mCPM, 'k', LineWidth=2); hold on; 
set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)

legend({'CS++' 'CS--' 'CS+-'})


%% 
md2p = mean(avCorrMM);
std2p = std(avCorrMM, [], 1); 
set2p = std2p ./ sqrt(size(avCorrMM, 1));
[h p ci ts] = ttest(avCorrMM);
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 



%plot(md2p); hold on; 
%plot(h)



times = 1:size(avCorrMM,2); 
shadedErrorBar(times, md2p, set2p, 'r', 1); hold on; 
plot (times, hb, 'Linewidth', 4)
%set(gca, 'xtick', [1 90 240], 'xticklabels', {'-1' '0' '1.5'}, 'xlim', [51 251])
%set(gca, 'FontSize', 12);
%plot([90 90],get(gca,'ylim'), 'k','lineWidth',1, 'Color', [.5 .5 .5]);
%plot(get(gca,'xlim'), [0 0 ], 'k','lineWidth',1, 'Color', [.5 .5 .5]);














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