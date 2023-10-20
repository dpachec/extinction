%%
%% Temporal RSA IN LOOP 
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast
clear , clc

listF2sav = {
                'POW_AMY_C_39-54_1_0_50-1_1_SISV-DISV';
                'POW_AMY_C_39-54_1_0_50-1_1_SISVAE-DISVAE';
        };   

t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        EEG = ALLEEG{subji};
        
        
        if ~isempty(EEG)

            
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); 

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = Ev2; 

            if strcmp(cfg.tyRSA, 'TR')
                EEG = normalize_baseline_EXT(EEG, [2501:3000]); 
                %EEG = normalize_EXT(EEG);  %across trials
                EEG = downsample_EEG_EXT(EEG); 
                cfg.oneListTraces = permute(EEG.data(:, 251:550,:), [3 1 2]); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'POW')
                EEG = extract_power_EXT(EEG, 0.01); 
                %EEG = normalize_baseline_EXT(EEG, [251:300]); 
                EEG = normalize_EXT(EEG);  %across trials
                cfg.oneListPow = EEG.power(:, :, : ,251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'PHA')
                EEG = normalize_EXT(EEG);  %across trials
                phaTS = extract_pha_EXT(EEG, cfg);
                cfg.oneListTraces = phaTS(:, :, 251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT3(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'PLV')
                EEG = normalize_EXT(EEG);  %across trials
                phaTS = extract_pha_EXT(EEG, cfg);
                cfg.oneListTraces = phaTS(:, :, 251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT5(out_contrasts, cfg);
                toc
            end
        
            ids{subji} = out_contrasts.allIDs; nnans{subji} = EEG.nan; 
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids', 'nnans');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%% plot 2 lines from TG 
clear
paths = load_paths_EXT; 
f2sav =   'POW_AMY_C_39-54_1_0_50-1_1_SICSPE-SICSME';
load ([ paths.results.rsa f2sav '.mat']);

ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

diff = cond1-cond2; 

for subji = 1:size(cond1, 1)
   cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
   cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
end       
cond1 = cond1B; cond2 = cond2B; 


d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1))
se2 = d2pstd2/sqrt(size(cond1, 1))

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
hb = h; hb(h==0) = nan; hb(hb==1) = -.03; 

times = (-.5:.01:2) + .25;
%times = 1:21
figure(); 
colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot(times, hb, LineWidth=6)
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 1);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 1);
%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%'ylim', [-.032 .035], 
title(f2sav, 'Interpreter','none')
exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)

%% permutations 2D (line plot)

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    %c1B = squeeze(cond1(:, 51:220)); 
    %c2B = squeeze(cond2(:, 51:220));
    c1B = squeeze(cond1(:, 51:151)); 
    c2B = squeeze(cond2(:, 51:151));
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :);
           c1B(subji, :) = c2B(subji, :);
           c2B(subji, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1B, c2B); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = min(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

clear p ratings2u mcsP
ratings2u = tObs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot TG
clear, clc
paths = load_paths_EXT; 

f2sav = 'POW_AMY_C_39-54_1_0_50-1_1_SISVA-DISVA';
%f2sav = 'TR_OFC_C_nan_0_0_50-1_1_SICSPE-SICSME';
                    

load ([ paths.results.rsa f2sav '.mat']);

ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
diff = cond1-cond2; 

[cond1 cond2] = rem_half_matrix(cond1, cond2);

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
%h(clustinfo.PixelIdxList{19}) = 1;

clim = [-.03 .03];
tRes = strsplit(f2sav, '_'); tRes = strsplit(tRes{7}, '-'); tRes = double(string(tRes{2}));
plot_TG_map(m1, m2, h, t, tRes, f2sav, clim)
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)




%% take mean in cluster 

for subji = 1:32

    c1 = squeeze(cond1(subji, :,:));
    c2 = squeeze(cond2(subji, :,:));

    mc1(subji,:) = mean(c1(clustinfo.PixelIdxList{4}));
    mc2(subji,:) = mean(c2(clustinfo.PixelIdxList{4}));

end



%% plot one bar
clc
ylim = [-0.1 0.1];
xlim = [0 3];
 
data.data = [mc1 mc2]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1 2], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'   '},     'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

%% plot 3 bars
clear, clc
load AMY_CSPE-CSMPP-CSMPM
ylim = [-0.1 0.1];
xlim = [0 4];
 
data.data = [mcCSPE mcCSMPP mcCSMPM]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1:3], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 3],'XTickLabel',{'1' '2' '3'},'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

%%
clear d4ANOVA
d4ANOVA = [mcCSPE ; mcCSMPP ; mcCSMPM];
d4ANOVA(:,2) = [ones(1,32) ones(1,32)*2 ones(1,32)*3];
d4ANOVA(any(isnan(d4ANOVA), 2), :) = [];
d4ANOVA(:,3) = [1:32 1:32 1:32];

x = RMAOV1(d4ANOVA);
[h p ci ts] = ttest(mcCSPE, mcCSMPP)
%[h p ci ts] = ttest(mcCSPE, mcCSMPM)

%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

%junts = cat(1, cond1(:, 3:21, 3:21), cond2(:, 3:21, 3:21));
%junts = cat(1, cond1(:, 4:18, 4:18), cond2(:, 4:18, 4:18));
%junts = cat(1, cond1(:, 4:13, 4:13), cond2(:, 4:13, 4:13));
%junts = cat(1, cond1(:, 51:end, 51:end), cond2(:, 51:end, 51:end));
%junts = cat(1, cond1(:, 51:151, 51:151), cond2(:, 51:151, 51:151));
%junts = cat(1, cond1(:, 26:225, 26:225), cond2(:, 26:225, 26:225));
junts = cat(1, cond1(:, 26:125, 26:125), cond2(:, 26:125, 26:125));

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
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
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

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)


%% PLV 200
clear
paths = load_paths_EXT; 
f2sav =   'PLV_OCC_V_3-8_0_0_200-10_1_SCA-DCA';

load ([ paths.results.rsa f2sav '.mat']);

ids = rem_nan_subj_EXT(out_rsa); 

cond1 = squeeze(out_rsa(:, 1, 6, 6)); 
cond2 = squeeze(out_rsa(:, 2, 6, 6)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

diff = cond1-cond2; 

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);


allC = [cond1 cond2];
boxplot(allC)

disp (['t: ' num2str(t) ' //  p = ' num2str(p)])

%% CHECK CONTEXT DURING ACQ AND EXT > GENERALIZATION
clear, clc
paths = load_paths_EXT; 
f2sav =   'POW_AMY_C_39-54_1_0_50-10_1_SCA-DCA';
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 
f2sav =   'POW_AMY_C_39-54_1_0_50-10_1_SCE-DCE';
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 




% % % % %  remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)
    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end
end

cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); cond1A(ids, :, :) = []; 
cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); cond2A(ids, :, :) = []; 
cond1E = squeeze(out_rsa_EXT(:, 1, :, :)); cond1E(ids, :, :) = []; 
cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); cond2E(ids, :, :) = []; 
diffA = cond1A-cond2A; 
diffE = cond1E-cond2E; 
m1 = squeeze(mean(diffA, 'omitnan')); 
m2 = squeeze(mean(diffE, 'omitnan')); 

[h p ci ts] = ttest(diffA, diffE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;

tRes = strsplit(f2sav, '_'); tRes = strsplit(tRes{7}, '-'); tRes = double(string(tRes{2}));
plot_TG_map(m1, m2, h, t, tRes, f2sav, [-.02 .02]); 
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)



%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(diffA, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, diffA(:, 4:23, 4:23), diffE(:, 4:23, 4:23));

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
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
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

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)


%% COMPUTE SUCCESIVE TRIALS SIMILARITY 
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast
clear , clc

listF2sav = {
                'POW_AMY_C_39-54_1_0_50-1_1_ALLE-TR';
                'POW_OFC_C_39-54_1_0_50-1_1_ALLE-TR';
                'POW_HPC_C_39-54_1_0_50-1_1_ALLE-TR';
                                
        };   

t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        EEG = ALLEEG{subji};
        
        
        if ~isempty(EEG)

            
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); %can't remove nan trials here or the order is lost

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = Ev2; 

            if strcmp(cfg.tyRSA, 'TR')
                EEG = normalize_baseline_EXT(EEG, [2501:3000]); 
                %EEG = normalize_EXT(EEG);  %across trials
                EEG = downsample_EEG_EXT(EEG); 
                cfg.oneListTraces = permute(EEG.data(:, 251:550,:), [3 1 2]); 
                out_contrasts = create_contrasts_trials_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'POW')
                EEG = extract_power_EXT(EEG, 0.01); 
                %EEG = normalize_baseline_EXT(EEG, [251:300]); 
                EEG = normalize_EXT(EEG);  %across trials
                cfg.oneListPow = EEG.power(:, :, : ,251:550); 
                out_contrasts = create_contrasts_trials_EXT(cfg);
                tic
                out_rsa{subji} = rsa_EXT6(out_contrasts, cfg);
                toc
            end
        
            ids{subji,:} = out_contrasts.allIDs;
            nnans{subji} = EEG.nan; 
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%% COMPUTE SUCCESIVE TRIALS SIMILARITY AND RATINGS CORRELATIONS
clear, clc
paths = load_paths_EXT; 

f2sav =  'POW_AMY_C_39-54_1_0_50-1_1_CSMPM-STR';
%f2sav =  'POW_HPC_C_39-54_1_0_50-1_1_ALLE-STR';

load ([ paths.results.rsa f2sav '.mat']);

load clustinfoRSA3

out_rsa = out_rsa(~cellfun('isempty', out_rsa));
ids = ids(~cellfun('isempty', ids));

for subji = 1:length(out_rsa)


    idsH = ids{subji}{1};
    %ratings = idsH(:, 7); 
    %ratings = idsH(:, 17) ;
    ratings = abs ( idsH(:, 7) - idsH(:,17) );
    %ratings = mean([idsH(:, 7) idsH(:,17)], 2); 
    

    orS = out_rsa{subji}; 
    %orS2 = squeeze(mean(mean(orS(:,:,3:23, 3:23), 4), 3));
    orS2 = squeeze(mean(mean(orS(:,:,26:100, 26:100), 4), 3))';
    
% %     clear orS2
% %     for triali = 1:size(orS, 2)
% %         orS2a = squeeze(orS(:,triali,:,:));
% %         orS2(triali,:) = mean(orS2a(clustinfo.PixelIdxList{11}));
% % 
% %         %dv = diag(orS2a); dv = dv(51:151);
% %         %orS2(triali,:) = mean(dv);
% %     end
% %     %allORS2(subji, :) = orS2;

    idN = isnan(orS2); 
    orS2(idN) = []; 
    ratings(idN) = []; 
    allRS(subji, :) = corr(orS2, ratings, 'type', 'k' ); 
        
% %         figure()
% %         %plot(ratings, orS2); hold on; 
% %         scatter(orS2, ratings); hold on; 
% %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% %         allSlopes(subji, :) = C(2);
% %         allIntercepts(subji, :) = C(1);
% %         %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
        

% %     figure()
% %     plot(orS2); hold on; 
% %     scatter(1:length(orS2), orS2); hold on; 
% %     h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% %     C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% %     allSlopes(subji, :) = C(2);
% %     allIntercepts(subji, :) = C(1);
% %     %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
% %     close all




end

%boxplot(allSlopes)
%[h p ci ts] = ttest(allSlopes); 

figure()
boxplot(allRS)
[h p ci ts] = ttest(allRS); 
disp(['T > ' num2str(ts.tstat) '  P > ' num2str(p)])



%% COMPUTE SIMILARITY OF EACH TRIAL TO ALL OTHER TRIALS
clear , clc

listF2sav = {
                'POW_OFC_C_39-54_1_0_50-1_1_ALLE-ATR';
                'POW_HPC_C_39-54_1_0_50-1_1_ALLE-ATR';
                'POW_AMY_C_39-54_1_0_50-1_1_SICSPE-ATR';
                'POW_AMY_C_39-54_1_0_50-1_1_SICSME-ATR';
                'POW_AMY_C_39-54_1_0_50-1_1_SICSMPP-ATR';
                'POW_AMY_C_39-54_1_0_50-1_1_SICSMPM-ATR';
                                
        };   

t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        EEG = ALLEEG{subji};
        
        
        if ~isempty(EEG)

            
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); %can't remove nan trials here or the order is lost

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = Ev2; 

            if strcmp(cfg.tyRSA, 'TR')
                EEG = normalize_baseline_EXT(EEG, [2501:3000]); 
                %EEG = normalize_EXT(EEG);  %across trials
                EEG = downsample_EEG_EXT(EEG); 
                cfg.oneListTraces = permute(EEG.data(:, 251:550,:), [3 1 2]); 
                out_contrasts = create_contrasts_trials_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'POW')
                EEG = extract_power_EXT(EEG, 0.01); 
                %EEG = normalize_baseline_EXT(EEG, [251:300]); 
                EEG = normalize_EXT(EEG);  %across trials
                cfg.oneListPow = EEG.power(:, :, : ,251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                out_rsa_p = rsa_EXT7(out_contrasts, cfg);
                id2u = out_contrasts.allIDs{1}(:,1);
                [id2unique idd2] = unique(id2u);
                
                clear allTRD
                for triali = 1:length(id2unique)
                    t2u = id2unique(triali);
                    m2c = id2u == t2u; 
                    allTRD(triali, :,:) = mean(out_rsa_p{1}(m2c, :, :));
                end
                out_rsa{subji,:} = allTRD; 
            end
        
            id2new{subji,:} = out_contrasts.allIDs{1}(idd2,:);
            nnans{subji,:} = EEG.nan; 
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'id2new');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end





%% COMPUTE SIMILARITY to all trials and AND RATINGS CORRELATIONS
clear, clc
paths = load_paths_EXT; 

f2sav =  'POW_AMY_C_39-54_1_0_50-1_1_ALLE-ATR';


load ([ paths.results.rsa f2sav '.mat']);

load clustinfoRSA3

out_rsa = out_rsa(~cellfun('isempty', out_rsa));
id2new = id2new(~cellfun('isempty', id2new));

for subji = 1:length(out_rsa)


    idsH = id2new{subji};
    ratings = idsH(:, 7); 
    
    orS = out_rsa{subji}; 
    %orS2 = squeeze(mean(mean(orS(:,:,3:23, 3:23), 4), 3));
    %orS2 = squeeze(mean(mean(orS(:,:,26:100, 26:100), 4), 3))';
    
    clear orS2
    for triali = 1:size(orS, 1)
        orS2a = squeeze(orS(triali,:,:));
        orS2(triali,:) = mean(orS2a(clustinfo.PixelIdxList{11}));

        %orS2a = squeeze(orS(triali,:,:));
        %dv = diag(orS2a); dv = dv(26:125);
        %orS2(triali,:) = mean(dv);
    end
    %allORS2(subji, :) = orS2;

    idM = isnan(ratings);
    idN = isnan(orS2); 
    orS2(idN | idM) = []; 
    ratings(idN | idM) = []; 
    allRS(subji, :) = corr(orS2, ratings, 'type', 'k' ); 
        
% %         figure()
% %         %plot(ratings, orS2); hold on; 
% %         scatter(orS2, ratings); hold on; 
% %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% %         allSlopes(subji, :) = C(2);
% %         allIntercepts(subji, :) = C(1);
% %         %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
        

% %     figure()
% %     plot(orS2); hold on; 
% %     scatter(1:length(orS2), orS2); hold on; 
% %     h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% %     C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% %     allSlopes(subji, :) = C(2);
% %     allIntercepts(subji, :) = C(1);
% %     %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
% %     close all




end

%boxplot(allSlopes)
%[h p ci ts] = ttest(allSlopes); 

figure()
boxplot(allRS)
[h p ci ts] = ttest(allRS); 
disp(['T > ' num2str(ts.tstat) '  P > ' num2str(p)])


%%
d2p = mean(allORS2, 'omitnan'); 

figure()
plot(d2p); hold on; 
scatter(1:length(d2p), d2p);
h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);




%%

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


clim = [-.03 .03];
tRes = strsplit(f2sav, '_'); tRes = strsplit(tRes{7}, '-'); tRes = double(string(tRes{2}));
plot_TG_map(m1, m2, h, t, tRes, f2sav, clim)
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)




%% take mean in cluster 

for subji = 1:32

    c1 = squeeze(cond1(subji, :,:));
    c2 = squeeze(cond2(subji, :,:));

    mc1(subji,:) = mean(c1(clustinfo.PixelIdxList{3}));
    mc2(subji,:) = mean(c2(clustinfo.PixelIdxList{3}));

end



%% plot one bar
clc
ylim = [-0.1 0.1];
xlim = [0 3];
 
data.data = [mc1 mc2]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1 2], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'   '},     'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

 
