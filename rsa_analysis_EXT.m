%% RSA ANALYSIS
%% Load data

clear 
paths = load_paths_EXT; 
%file2load = ['allS_' 'Amygdala' '_C'];
file2load = ['allS_' 'inferiortemporal_middletemporal_superiortemporal_bankssts_fusiform_temporalpole_lateraloccipital_lingual_parahippocampal_cuneus_pericalcarine' '_C'];
load ([paths.results.power file2load]); 


%%contrast based RSA

%freqs_avTimeFeatVect_freqResolv(0-1)_trials/noTrials_win-width_mf
clc
clearvars -except ALLEEG paths file2load

%f2sav = '3-8_1_0_0_50-1_DISVA-DIDVA_TG'; 
f2sav = '39-54_1_0_0_50-1_1_SICSPA-SICSMA'; 
%f2sav = '3-54_1_0_0_50-1_1_DISCA-DIDCA-SICSPE-SICSME-DISVA-DIDVA_TG'; 

cfg = getParams_EXT(f2sav);

t1 = datetime; 
for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        
        cfg.oneListIds = Ev2; 
        cfg.oneListPow = EEG.power(:, :, : ,251:470); 

        out_contrasts = create_contrasts_EXT(cfg);

        out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
        
        
        
    end

end


mkdir ([paths.results.rsa]);
save([ paths.results.rsa f2sav '_' file2load '.mat'], 'out_rsa');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% 
clearvars -except ALLEEG f2sav paths file2load
f2sav = [ '39-54_1_0_0_50-1_1_SICSPE-SICSME_' file2load ]; 
load([paths.results.rsa f2sav '.mat']);



%% remove hack 
ids = []; 
for subji = 1:length(ALLEEG)

    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end

end

%%
cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

[max2u id] = max(abs(allSTs));
tObs = allSTs(id); 


tiledlayout(1,3);
nexttile
%imagesc(m1);  axis square
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; %colorbar
nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square; 
nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; %colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);



%axesHandles = findall(0, 'type', 'axis');
%set(axesHandles,'equal'); 



%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

%junts = [avSICP; avSICM];
junts = cat(1, cond1, cond2);

clear max_clust_sum_perm
for permi = 1:nPerm
    
    [M,N] = size(realCondMapping);
    rowIndex = repmat((1:M)',[1 N]);
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);

    cond1P = junts(fakeCondMapping == 0, :,:);
    cond2P = junts(fakeCondMapping == 1, :,:);

    diffC = cond1P - cond2P; 
    diffC = diffC(:,51:170, 51:170);
    [h p ci ts] = ttest(diffC); 
    t = ts.tstat; 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

%% 


%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)







%%



















%% count channels 

clear nChans
for subji = 1:length(ALLEEG1)
    EEG = ALLEEG1{subji};
    if~isempty(EEG)
        nChans(subji,:) = length(EEG.chanlocs);
    end
end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])

%% select first 2 in each hemisphere (Amygdala)  

clearvars -except ALLEEG paths file2load
clc 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 

        
        ids2rem1 = []; 
        %ids2rem1 = contains(chans, 'Hippocampus')
        %ids2rem1 = contains(chans, 'Right')
        %ids2rem = logical(ids2rem1+ids2rem2)
        ids2rem = ids2rem1; 
        EEG.chanlocs(ids2rem) = []; 
        if size(EEG.chanlocs, 2) >0 & size(EEG.chanlocs, 1) >0 
            EEG.power(:, ids2rem, :, :) = []; 
            ALLEEG1{subji,:} = EEG; 
        end
        
%         hy1 = find(contains(chans, 'Left')); 
%         if ~isempty(hy1) idL = hy1(1); end
%         hy2 = find(contains(chans, 'Right'));
%         if ~isempty(hy2) idR = hy2(1); end
%         
%         if ~isempty(hy1) idx = idL; end
%         if ~isempty(hy2) idx = idR; end
%         if ~isempty(hy1) & ~isempty(hy2) idx = [idL idR]; end
%         
%         EEG.power = EEG.power(:, idx, :, :); 
%         EEG.chanlocs = EEG.chanlocs(idx);
%         ALLEEG1{subji,:} = EEG; 



    end


    

end



%% compute neural RDM 
%freqs_avTimeFeatVect_freqResolv(0-1)_trials/noTrials_win-width_mf
clc
clearvars -except ALLEEG ALLEEG1 paths file2load

f2sav = '3-8_1_0_0_50-1_DISVA-DIDVA_TG'; 
cfg = getParams_EXT(f2sav);


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        ids1 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '1');
        ids2 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '2');
        ids3 =  strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        
        all1 = EEG.power(ids1, :, : ,201:550); 
        all2 = EEG.power(ids2, :, : ,201:550); 
        all3 = EEG.power(ids3, :, : ,201:550); 

        allIts = cat(1, all1, all2, all3);
        neuralRDMALL = createNeuralRDMs_EXT(cfg, allIts);
        
        neuralRDMALL = remDiagRDM_TS(neuralRDMALL); 

        %neuralRDMALL2p = mean(neuralRDMALL(:,:,145:155), 3); 
        %nanTr2 = find(all(isnan(neuralRDMALL2p), 2)); %rows id if there are nans 
        %allNTR(subji,:) = length(nanTr2); 


        % % % % plot RDMS
        %if length(nanTr2) < 24
            %neuralRDMAllSort(subji, :, :) = neuralRDMALL2p;
            %figure()
            %imagesc(neuralRDMALL2p)
        %end
        

        sRDM = size(neuralRDMALL, 1); 
        sRDM1= size(all1, 1);
        sRDM2= size(all2, 1);

            
    % % %             % % % CHECK THAT TEMPLATE USED ARE CORRECT 
    % % %             SI1_T = zeros(sRDM);SI2_T = zeros(sRDM);SI3_T = zeros(sRDM);
    % % %             SI1_T(1:sRDM1,1:sRDM1) = 1; 
    % % %             SI2_T(sRDM1+1:sRDM2+sRDM1,sRDM1+1:sRDM2+sRDM1) = 1; 
    % % %             SI3_T(sRDM1+sRDM2+1:end,sRDM1+sRDM2+1:end) = 1; 
    % % %             figure; imagesc(SI1_T); axis square
    % % %             figure; imagesc(SI2_T); axis square
    % % %             figure; imagesc(SI3_T); axis square
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


            %avCorrSI1TR = mean(mean(SI1(:, :,145:155), 1, 'omitnan'), 3, 'omitnan');
            %avCorrSI2TR = mean(mean(SI2(:, :,145:155), 1, 'omitnan'), 3, 'omitnan');
            %avCorrSI3TR = mean(mean(SI3(:, :,145:155), 1, 'omitnan'), 3, 'omitnan');
            %avCorrTRALL(subji, :)= [avCorrSI1TR avCorrSI2TR avCorrSI3TR];
            


        %end
       
    end

end

% % check if subjects have nans in one condition and remove from the others
cc1 = find(any(isnan(avCorrSI1), 2)   | avCorrSI1(:, 1)==0); %rows id if there are nans 
cc2 = find(any(isnan(avCorrSI2), 2)   | avCorrSI2(:, 1)==0);
cc3 = find(any(isnan(avCorrSI3), 2)   | avCorrSI3(:, 1)==0);
ccsv = find(any(isnan(avCorrDISV), 2) | avCorrDISV(:, 1)==0);
ccdv = find(any(isnan(avCorrDIDV), 2) | avCorrDIDV(:, 1)==0);

sub2rem = unique([cc1;cc2;cc3;ccsv;ccdv]); 

avCorrSI1(sub2rem, : ) = [];
avCorrSI2(sub2rem, : ) = [];
avCorrSI3(sub2rem, : ) = [];
avCorrDISV(sub2rem, : ) = [];
avCorrDIDV(sub2rem, : ) = [];


disp('done');







%% plot 2 conditions acquisition

sub2exc = [];

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

itSpec = avSICP - avDISV; 
maITSPEC = squeeze(mean(itSpec));
[h p ci ts] = ttest(itSpec); 
t = ts.tstat; 

times = (-1:.01:2) + .25
%times = 1:size(mavSICP,2); 
figure(); 
plot(times, mavSICP, 'r', LineWidth=2); hold on; 
plot(times, mavSICM, 'k', LineWidth=2); hold on; 
set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)
legend({'CS+' 'CS-'})

figure(); 
plot(times, mavDISV, 'r', LineWidth=2); hold on; 
plot(times, mavDIDV, 'k', LineWidth=2); hold on; 
set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)
legend({'DISV' 'DIDV'})

figure(); 
plot(times, maITSPEC, 'k', LineWidth=2); hold on; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.002; 
%plot (times, hb, LineWidth=5)
scatter(times, hb)
set(gca, 'xlim' ,[-.5 2], 'Fontsize', 18)
legend({'Item-specific activity'})



%% plot SI CS+ SI CS- conditions acquisition variance 

mSICP = mean(avSICP);
stdSICP = std(avSICP, [], 1); 
seSICP = stdSICP / sqrt(size(avSICP, 1));

mSICM = mean(avSICM);
stdSICM = std(avSICM, [], 1); 
seSICM = stdSICM / sqrt(size(avSICM, 1));

diffC = avSICP - avSICM; 
[h p ci ts] = ttest(diffC); 
t = ts.tstat; 
clustinfo = bwconncomp(h);

clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObsCSPM= allSTs(id) ;
end

%h = zeros(1, size(avSICP, 2));
h (1:100) = 0



hb = h; hb(h==0) = nan; hb(hb==1) = -.02; 


colors2use = brewermap([4],'*Set1')*0.75;

%times = 1:size(avSICP,2); 
times = (-1:.01:2) + .25;
shadedErrorBar(times, mSICP, seSICP, {'Color', colors2use(1, :)}, 1); hold on; 
shadedErrorBar(times, mSICM, seSICM, {'Color', colors2use(2, :)}, 1); hold on; 
xlabel('Time (s)')
ylabel('Rho')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)






%% plot 2 conditions EXTINCTION

sub2exc = [4 13];

avSI1 = avCorrSI1; 
avSICP = avCorrSI1; 
avSICM = squeeze(mean(cat(3, avCorrSI2, avCorrSI3), 3)); 
avDISV = avCorrDISV; avDIDV = avCorrDIDV;

avSICM(sub2exc, :) = []; avSICP(sub2exc, :) = []; avDISV(sub2exc, :) = []; avDIDV(sub2exc, :) = []; 

mavSICP = mean(avSICP); mavSICM = mean(avSICM); mavDISV = mean(avDISV); mavDIDV = mean(avDIDV); 

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
legend({'DISV' 'DIDV'})

%% plot SI CS+ SI CS- conditions EXTINCTION variance 

mSICP = mean(avSICP);
stdSICP = std(avSICP, [], 1); 
seSICP = stdSICP / sqrt(size(avSICP, 1));

mSICM = mean(avSICM);
stdSICM = std(avSICM, [], 1); 
seSICM = stdSICM / sqrt(size(avSICM, 1));

diffC = avSICP - avSICM; 
[h p ci ts] = ttest(diffC); 
t = ts.tstat; 
clustinfo = bwconncomp(h);

clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObsCSPM= allSTs(id) ;
end


%h = zeros(1, size(avSICP, 2));
h (200:end) = 0

%times = 1:size(avSICP,2); 
times = (-1:.01:2) + .25;
hb = h; hb(h==0) = nan; hb(hb==1) = -.02; 

colors2use = brewermap([4],'*Set1')*0.75;
shadedErrorBar(times, mSICP, seSICP, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, mSICM, seSICM,  {'Color',colors2use(2,:)}, 1); hold on; 

xlabel('Time (s)')
ylabel('Rho')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)


%% plot VALENCE WITH VARIANCE

mDISV = mean(avDISV);
stdDISV = std(avDISV, [], 1); 
seDISV = stdDISV / sqrt(size(avDISV, 1));

mDIDV = mean(avDIDV);
stdDIDV = std(avDIDV, [], 1); 
seDIDV = stdDIDV / sqrt(size(avDIDV, 1));

diffC = avDISV - avDIDV; 
[h p ci ts] = ttest(diffC); 
t = ts.tstat; 
clustinfo = bwconncomp(h);

clear allSTs tObs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(allSTs);
    tObs= allSTs(id) 
    tObsVAL = tObs; 
end

hb = h; hb(h==0) = nan; hb(hb==1) = -.01; 


times = (-1:.01:2) + .25;

colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  mDISV, seDISV, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, mDIDV, seDIDV,  {'Color',colors2use(2,:)}, 1); hold on; 

xlabel('Time (s)')
ylabel('Rho')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)
%set(gca, 'xtick', [1 90 240], 'xticklabels', {'-1' '0' '1.5'}, 'xlim', [51 251])
%plot([90 90],get(gca,'ylim'), 'k','lineWidth',1, 'Color', [.5 .5 .5]);
%plot(get(gca,'xlim'), [0 0 ], 'k','lineWidth',1, 'Color', [.5 .5 .5]);


%% plot 3 items separately

sub2exc = []

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


%% plot 3 items separately with variance

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



%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(avDISV, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

%junts = [avSICP; avSICM];
junts = [avDISV; avDIDV];

clear max_clust_sum_perm
for permi = 1:nPerm
    
    [M,N] = size(realCondMapping);
    rowIndex = repmat((1:M)',[1 N]);
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);


    cond1 = junts(fakeCondMapping == 0, :);
    cond2 = junts(fakeCondMapping == 1, :);

    diffC = cond1 - cond2; 
    [h p ci ts] = ttest(diffC); 
    t = ts.tstat; 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

%% 

%tObs = tObsCSPM;
tObs = tObsVAL; 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)










%% CONTEXT ANALYSIS
%freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf
clc
clearvars -except ALLEEG ALLEEG1 paths file2load

f2sav = '39-54_1_0_0_50-1'; 
cfg = getParams_EXT(f2sav);


for subji = 1:length(ALLEEG1)
    
    EEG = ALLEEG1{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        ctxs = Ev2(:, 3); 
        ctxs = cellfun(@(x) x(2), ctxs, 'un', 0);

        ids1 =  strcmp(Ev2(:, 2), '1') & strcmp(ctxs, '1');
        ids2 =  strcmp(Ev2(:, 2), '1') & strcmp(ctxs, '2');
        ids3 =  strcmp(Ev2(:, 2), '1') & strcmp(ctxs, '3');
        ids4 =  strcmp(Ev2(:, 2), '1') & strcmp(ctxs, '4');

        
        all1 = EEG.power(ids1, :, : ,201:550); 
        all2 = EEG.power(ids2, :, : ,201:550); 
        all3 = EEG.power(ids3, :, : ,201:550); 
        all4 = EEG.power(ids4, :, : ,201:550); 

        allIts = cat(1, all1, all2, all3, all4);
        neuralRDMALL = createNeuralRDMs_EXT(cfg, allIts);
        neuralRDMALL = remDiagRDM_TS(neuralRDMALL); 


        sRDM = size(neuralRDMALL, 1); 
        sRDM1= size(all1, 1);
        sRDM2= size(all2, 1);
        sRDM3 = size(all3, 1);
        sRDM4 = size(all4, 1);

        SI1_T = zeros(sRDM);SI2_T = zeros(sRDM);SI3_T = zeros(sRDM);SI4_T = zeros(sRDM);S_ALL = zeros(sRDM);
        SI1_T(1:sRDM1,1:sRDM1) = 1; 
        SI2_T(sRDM1+1:sRDM2+sRDM1,sRDM1+1:sRDM2+sRDM1) = 1; 
        SI3_T(sRDM1+sRDM2+1:sRDM1+sRDM2+sRDM3,sRDM1+sRDM2+1:sRDM1+sRDM2+sRDM3) = 1; 
        SI4_T(sRDM1+sRDM2+sRDM3+1:end,sRDM1+sRDM2+sRDM3+1:end) = 1; 

% % % %         % % % CHECK THAT TEMPLATE USED ARE CORRECT 
% % % %         figure; imagesc(SI1_T); axis square
% % % %         figure; imagesc(SI2_T); axis square
% % % %         figure; imagesc(SI3_T); axis square
% % % %         figure; imagesc(SI4_T); axis square
% % % %         figure; imagesc(S_ALL); axis square

    
        S_ALL(SI1_T == 1) = 1;S_ALL(SI2_T == 1) = 1;S_ALL(SI3_T == 1) = 1;S_ALL(SI4_T == 1) = 1;
       
        parfor timei = 1:size(neuralRDMALL, 3)
            tmpM = squeeze(neuralRDMALL(:, :, timei)); 
            SC(timei) = mean(tmpM(S_ALL == 1), 'omitnan');
            DC(timei) = mean(tmpM(S_ALL == 0), 'omitnan');
        end
        avSC(subji, :) = SC;
        avDC(subji, :) = DC;
        
    end

end

% % check if subjects have nans in one condition and remove from the others
cc1 = find(any(isnan(avSC), 2)); %rows id if there are nans 
cc2 = find(any(isnan(avDC), 2)); %rows id if there are nans 

sub2rem = [cc1 cc2]; 
avSC(sub2rem, : ) = [];
avDC(sub2rem, : ) = [];


avSC =  avSC (any(avSC ,2),:);
avDC =  avDC (any(avDC ,2),:);

disp('done');



%% plot SAME CONTEXT vs DIFF CONTEXTs

mSC = mean(avSC);
stdSC = std(avSC, [], 1); 
seSC = stdSC / sqrt(size(avSC, 1));

mDC = mean(avDC);
stdDC = std(avDC, [], 1); 
seDC = stdDC / sqrt(size(avDC, 1));

diffC = avSC - avDC; 
[h p ci ts] = ttest(diffC); 
t = ts.tstat; 
clustinfo = bwconncomp(h);

clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObsCSPM= allSTs(id) ;
end


%h = zeros(1, size(avSC, 2));
%h (200:end) = 0

%times = 1:size(avSICP,2); 
times = (-1:.01:2) + .25;
hb = h; hb(h==0) = nan; hb(hb==1) = -.02; 
colors2use = brewermap([4],'*Set1')*0.75;
shadedErrorBar(times, mSC, seSC,{'Color',colors2use(4,:)}, 1); hold on; 
shadedErrorBar(times, mDC, seDC,{'Color',colors2use(3,:)}, 1); hold on; 


xlabel('Time (s)')
ylabel('Rho')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);


exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)




%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(avSC, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = [avSC; avDC];
%junts = [avDISV; avDIDV];

clear max_clust_sum_perm
for permi = 1:nPerm
    
    [M,N] = size(realCondMapping);
    rowIndex = repmat((1:M)',[1 N]);
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);


    cond1 = junts(fakeCondMapping == 0, :);
    cond2 = junts(fakeCondMapping == 1, :);

    diffC = cond1 - cond2; 
    [h p ci ts] = ttest(diffC); 
    t = ts.tstat; 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

%% 
%tObs =  -30.4546%-86.4470;
tObs = tObsCSPM;
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)







%% plot examples 
subj = 13
d2p = squeeze(neuralRDMAllSort(subj, :, :));
figure()
imagesc(d2p); axis square

%% plot GA RDM

figure()
d2p = squeeze(mean(neuralRDMAllSort, 1, 'omitnan'))
imagesc(d2p); axis square; colorbar






%%