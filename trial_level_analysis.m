%% Correlate amygdala power IN CLUSTER with different trial-level metrics
%% 
clear, clc

%tp2use = 36:46;
%tp2use = 20:32; % check OFC
tp2use = 33:41; % > AMY theta power effect


paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_1-44_1_0_500-50'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_1-44_1_0_500-50'])



minTr = 8; 


if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end

sub2exc = [37];

allRHO = nan(50, 1); 

for subji = 1:50

    amyPOW = allPOWAMY{subji, 1}; 
    amyPOWIDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 

    if ~isempty(amyPOW) & ~isempty(rsa2T)
       
        [C i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW(i1, :); 
        amyPOWIDs = amyPOWIDs(i1,:); 
        rsa2T = mean(rsa2T(i2, tp2use), 2); 


        % % % z-score amygdala power separately for CS+ and CS-
        amyPOWCSp = amyPOW(amyPOWIDs(:, 8) == 1); 
        amyPOWCSm = amyPOW(amyPOWIDs(:, 8) == 0); 
        amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
        amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
        amyPOW(amyPOWIDs(:, 8) == 1) = amyPOWCSp; 
        amyPOW(amyPOWIDs(:, 8) == 0) = amyPOWCSm; 

        % % % z-score RSA2T metric separately for CS+ and CS-
        % rsa2TCSp = rsa2T(rsa2TIDs(:, 8) == 1); 
        % rsa2TCSm = rsa2T(rsa2TIDs(:, 8) == 0); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2T(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        % rsa2T(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nanIds = isnan(amyPOW); 
        amyPOW(nanIds) = []; 
        rsa2T(nanIds) = []; 


        % [B, tF1] = rmoutliers(rsa2T, 'percentiles', [10 90]); 
        % [B, tF2] = rmoutliers(amyPOW, 'percentiles', [10 90]); 
        % [B, tF1] = rmoutliers(rsa2T, 'mean'); 
        % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
        % tF3 = tF1 | tF2; 
        % nOut(subji, :) = sum(tF3); 
        % rsa2T(tF3) = []; 
        % amyPOW(tF3) = []; 

        % % % 
        % % % figure()
        % % % scatter(amyPOW, rsa2T, 150, 'filled');
        % % % h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
        % % % C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
        % % % allSlopes(subji, :) = C(2);
        % % % allIntercepts(subji, :) = C(1);
        % % % set(gca, 'ylim', [-.15 .15], 'xlim', [-3 3], 'Fontsize', 24)
        % % % %set(gca, 'ylim', [-3 3], 'xlim', [-3 3], 'Fontsize', 24)


        if length(amyPOW) > minTr % at least five trials for the correlation analysis
            allRHO(subji, :) = corr(amyPOW, rsa2T, 'type', 's');
        end

    else
        
        allRHO(subji, :) = nan; 

    end

end

%allRHO(sub2exc) = nan; 
allRHO(isnan(allRHO)) = []; 
[h p ci ts] = ttest(atanh(allRHO));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


%% plot one bar
clear data
data.data = atanh(allRHO); 

figure(); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 35, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
%set(gca, 'ylim', [-.75 .75])
set(gca, 'ylim', [-1 .2])
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'myP.png', 'Resolution', 300);


%% Correlate amygdala power IN CLUSTER with different trial-level metrics FOR DIFFERENT CONDITIONS
clearvars -except allRHO_1 allRHO_2 allRHO_3
clc

%tp2use = 36:46;
%tp2use = 20:32; % check OFC

con2u = 3; 
tp2use = 33:41; % > AMY theta power effect


paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_9-54_1_0_500-100'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OCC_CE_3-54_1_0_500-100'])


minTr = 8; 


if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end

sub2exc = [37];

allRHO = nan(50, 1); 

for subji = 1:50

    amyPOW = allPOWAMY{subji, 1}; 
    amyPOWIDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 

    if ~isempty(amyPOW) & ~isempty(rsa2T)
        
        conids = rsa2TIDs(:, 6) == con2u; 
        amyPOW = amyPOW(conids); 
        amyPOWIDs = amyPOWIDs(conids,:); 
        
        [C i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW(i1, :); 
        amyPOWIDs = amyPOWIDs(i1,:); 
        rsa2T = mean(rsa2T(i2, tp2use), 2); 



        % % % z-score amygdala power separately for CS+ and CS-
        amyPOWCSp = amyPOW(amyPOWIDs(:, 8) == 1); 
        amyPOWCSm = amyPOW(amyPOWIDs(:, 8) == 0); 
        amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
        amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
        amyPOW(amyPOWIDs(:, 8) == 1) = amyPOWCSp; 
        amyPOW(amyPOWIDs(:, 8) == 0) = amyPOWCSm; 

        % % % z-score RSA2T metric separately for CS+ and CS-
        % % rsa2TCSp = rsa2T(rsa2TIDs(:, 8) == 1); 
        % % rsa2TCSm = rsa2T(rsa2TIDs(:, 8) == 0); 
        % % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % % rsa2T(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        % % rsa2T(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nanIds = isnan(amyPOW); 
        amyPOW(nanIds) = []; 
        rsa2T(nanIds) = []; 


        % [B, tF1] = rmoutliers(rsa2T, 'percentiles', [10 90]); 
        % [B, tF2] = rmoutliers(amyPOW, 'percentiles', [10 90]); 
        % [B, tF1] = rmoutliers(rsa2T, 'mean'); 
        % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
        % tF3 = tF1 | tF2; 
        % nOut(subji, :) = sum(tF3); 
        % rsa2T(tF3) = []; 
        % amyPOW(tF3) = []; 

        % % % 
        % % % figure()
        % % % scatter(amyPOW, rsa2T, 150, 'filled');
        % % % h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
        % % % C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
        % % % allSlopes(subji, :) = C(2);
        % % % allIntercepts(subji, :) = C(1);
        % % % set(gca, 'ylim', [-.15 .15], 'xlim', [-3 3], 'Fontsize', 24)
        % % % %set(gca, 'ylim', [-3 3], 'xlim', [-3 3], 'Fontsize', 24)


        if length(amyPOW) > minTr % at least five trials for the correlation analysis
            allRHO(subji, :) = corr(amyPOW, rsa2T, 'type', 's');
        end

    else
        
        allRHO(subji, :) = nan; 

    end

end

%%
clc
%allRHO(sub2exc) = nan; 
%allRHO(isnan(allRHO)) = []; 
[h p ci ts] = ttest(atanh(allRHO_1), atanh(allRHO_3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


%%
%allRHO(sub2exc) = nan; 
allRHO(isnan(allRHO)) = []; 
[h p ci ts] = ttest(atanh(allRHO));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

%%
clc 
isnan1 = isnan(allRHO_1); 
isnan2 = isnan(allRHO_2); 
isnan3 = isnan(allRHO_3); 
allisN = isnan1|isnan2|isnan3; 

allRHO_1(allisN) = []; 
allRHO_2(allisN) = []; 
allRHO_3(allisN) = []; 

data = [allRHO_1, allRHO_2, allRHO_3];
nSubj = size(allRHO_1, 1); 
d4anova = data(:);
d4anova = d4anova(:); 
d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];

[p f] = RMAOV1(d4anova);

boxplot(data)

%%

[h p ci ts] = ttest(allRHO_1, allRHO_3)


%% plot one bar
clear data
data.data = atanh(allRHO); 

figure(); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 35, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
%set(gca, 'ylim', [-.75 .75])
set(gca, 'ylim', [-1 .2])
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'myP.png', 'Resolution', 300);



%% Correlate amygdala POWER IN CLUSTER with different trial-level metrics (ALL TIME POINTS)
clear, clc

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_1-44_1_0_500-50'])



load ([paths.results.trial_based 'trlSTA_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_1-44_1_0_500-50'])

sub2exc = [37]; 
minTr = 5; 



if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end



win_width = 10; 
mf = 1; 
nTimepoints = 51; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins); 

for subji = 1:50

    amyPOW1 = allPOWAMY{subji, 1}; 
    amyPOWIDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 
    
    if ~isempty(amyPOW1) & ~isempty(rsa2T)

        [C i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW1(i1, :); 
        amyPOWIDsh = amyPOWIDs(i1,:); 
        rsa2T = rsa2T(i2, :); 

        nanIds = isnan(amyPOW); 
        amyPOW(nanIds,:) = []; 
        rsa2T(nanIds,:) = []; 
        amyPOWIDsh(nanIds,:) = []; 
        rsa2TIDs(nanIds,:) = []; 


        % % % z-score amygdala power separately for CS+ and CS-
        amyPOWCSp = amyPOW(amyPOWIDsh(:, 8) == 1); 
        amyPOWCSm = amyPOW(amyPOWIDsh(:, 8) == 0); 
        amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
        amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
        amyPOW(amyPOWIDsh(:, 8) == 1) = amyPOWCSp; 
        amyPOW(amyPOWIDsh(:, 8) == 0) = amyPOWCSm; 
        % % % z-score RSA2T metric separately for CS+ and CS-
        % rsa2TCSp = rsa2TT(rsa2TIDs(:, 8) == 1); 
        % rsa2TCSm = rsa2TT(rsa2TIDs(:, 8) == 0); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2TT(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        % rsa2TT(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nTrials (subji, : ) = length(amyPOW); 
        



         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT = mean(rsa2T(:, timeBins), 2);     
    
            % [B, tF1] = rmoutliers(rsa2TT, 'percentiles', [5 95]); 
            % [B, tF2] = rmoutliers(amyPOW, 'percentiles', [5 95]); 
            % [B, tF1] = rmoutliers(rsa2TT, 'mean'); 
            % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
            % tF3 = tF1 | tF2; 
            % nOut(subji, timei) = sum(tF3); 
            % rsa2TT(tF3,:) = []; 
            % amyPOWH = amyPOW; 
            % amyPOWH(tF3,:) = []; 

            if length(amyPOW) > minTr % at least five trials for the correlation analysis
                allRHO(subji, timei, :) = corr(amyPOW, rsa2TT, 'type', 's');
            
            else
            
                allRHO(subji, timei, :) = nan; 
            end
    
         end

    end
end

allRHO(sub2exc, :) = []; 

allRHO = allRHO(any(allRHO,2),:);
[h p ci ts] = ttest(atanh(allRHO));
t = squeeze(ts.tstat); 

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end


%h(1:10) = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

mAR = mean(allRHO);
stdAR = std(allRHO); 
seAR = stdAR / sqrt(size(allRHO, 1))

xStart = -.2; dx = 0.05; times = xStart + (0:bins-1)*dx;
h(2, :) = times; 
figure()
%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'Set3')*0.75;
shadedErrorBar(times, mAR, seAR, {'Color',colors2use(1,:)}, 1); hold on; 

set(gca, 'ylim', [-.4 .2], 'xlim', [-.25 1.75]) % 

plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(times, hb, Linewidth=6); 

set(gca, 'Fontsize', 30)

exportgraphics(gcf, ['myP.png'], 'Resolution',150)







%% Permutations
clearvars -except max_clust_obs allPOWAMY cond2u
clc


nPerm = 1000; 
t4P = 3:20; 
nTimepoints = 18; %%all2all file is stored within this cell array
win_width = 5; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

for permi = 1:nPerm

    clearvars -except permi allPOWAMY cond2u bins mf win_width allallRHO max_clust_perm max_clust_obs nPerm t4P

    for subji = 1:50

        amyPOW1 = allPOWAMY{subji, 1}; 
        amyPOWIDs = double(string(allPOWAMY{subji, 2})); 
    
        rsa2T = cond2u{subji, 1}; 
        rsa2TIDs = cond2u{subji, 2}; 
        
        if ~isempty(amyPOW1) & ~isempty(rsa2T)
            [i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
            amyPOW = amyPOW1(i2, :); 
            amyPOWIDsh = amyPOWIDs(i2,:); 
        
            nanIds = isnan(amyPOW); 
            amyPOW(nanIds,:) = []; 
            rsa2T(nanIds,:) = []; 
            amyPOWIDsh(nanIds,:) = []; 
            rsa2TIDs(nanIds,:) = []; 


            ids4perm = randperm(size(rsa2TIDs, 1));
            rsa2T = rsa2T(ids4perm, t4P); 
            rsa2TIDs = rsa2TIDs(ids4perm, :); 
            

    
             for timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    
                
    
                rsa2TT = mean(rsa2T(:, timeBins), 2); 
        
                % % % z-score amygdala power separately for CS+ and CS-
                amyPOWCSp = amyPOW(amyPOWIDsh(:, 8) == 1); 
                amyPOWCSm = amyPOW(amyPOWIDsh(:, 8) == 0); 
                amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
                amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
                amyPOW(amyPOWIDsh(:, 8) == 1) = amyPOWCSp; 
                amyPOW(amyPOWIDsh(:, 8) == 0) = amyPOWCSm; 
                % % % z-score RSA2T metric separately for CS+ and CS-
                % rsa2TCSp = rsa2TT(rsa2TIDs(:, 8) == 1); 
                % rsa2TCSm = rsa2TT(rsa2TIDs(:, 8) == 0); 
                % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
                % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
                % rsa2TT(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
                % rsa2TT(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 
                % 
        
                if length(amyPOW) > 5 % at least five trials for the correlation analysis
                    allRHO(subji, timei, :) = corr(amyPOW, rsa2TT, 'type', 's');
                else
                    allRHO(subji, timei, :) = nan; 
                end
        
            end
        end
    end

    allRHO = allRHO(any(allRHO,2),:);
    [h p ci ts] = ttest(atanh(allRHO));
    t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, :) = allSTs(id); 
    else
        max_clust_perm(permi, :) = 0; 
    end

    %allallRHO{permi,:} = allRHO; 

end





allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm
















%% Correlate amygdala POWER IN CLUSTER with different trial-level metrics SEPARATELY FOR EACH CONDITION
clear, clc


con2u = 1; 

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_1-44_1_0_500-50'])



load ([paths.results.trial_based 'trlSTA_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_1-44_1_0_500-50'])

sub2exc = [37]; 
minTr = 5; 



if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end



win_width = 10; 
mf = 1; 
nTimepoints = 51; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins); 

for subji = 1:50

    amyPOW1 = allPOWAMY{subji, 1}; 
    amyPOWIDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 
    
    if ~isempty(amyPOW1) & ~isempty(rsa2T)
        [i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW1(i2, :); 
        amyPOWIDsh = amyPOWIDs(i2,:); 
    
        nanIds = isnan(amyPOW); 
        amyPOW(nanIds,:) = []; 
        rsa2T(nanIds,:) = []; 
        amyPOWIDsh(nanIds,:) = []; 
        rsa2TIDs(nanIds,:) = []; 


        % % % z-score amygdala power separately for CS+ and CS-
        amyPOWCSp = amyPOW(amyPOWIDsh(:, 8) == 1); 
        amyPOWCSm = amyPOW(amyPOWIDsh(:, 8) == 0); 
        amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
        amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
        amyPOW(amyPOWIDsh(:, 8) == 1) = amyPOWCSp; 
        amyPOW(amyPOWIDsh(:, 8) == 0) = amyPOWCSm; 
        % % % z-score RSA2T metric separately for CS+ and CS-
        % rsa2TCSp = rsa2T(rsa2TIDs(:, 8) == 1); 
        % rsa2TCSm = rsa2T(rsa2TIDs(:, 8) == 0); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2T(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        % rsa2T(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nTrials (subji, : ) = length(amyPOW); 
        
        conids = rsa2TIDs(:, 6) == con2u; 
        amyPOW = amyPOW(conids); 
        rsa2T = rsa2T(conids, :); 
        
         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT = mean(rsa2T(:, timeBins), 2);     
    
            % [B, tF1] = rmoutliers(rsa2TT, 'percentiles', [5 95]); 
            % [B, tF2] = rmoutliers(amyPOW, 'percentiles', [5 95]); 
            % [B, tF1] = rmoutliers(rsa2TT, 'mean'); 
            % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
            % tF3 = tF1 | tF2; 
            % nOut(subji, timei) = sum(tF3); 
            % rsa2TT(tF3,:) = []; 
            % amyPOWH = amyPOW; 
            % amyPOWH(tF3,:) = []; 

            if length(amyPOW) > minTr % at least five trials for the correlation analysis
                allRHO(subji, timei, :) = corr(amyPOW, rsa2TT, 'type', 's');
            
            else
            
                allRHO(subji, timei, :) = nan; 
            end
    
         end

    end
end

allRHO(sub2exc, :) = []; 

allRHO = allRHO(any(allRHO,2),:);
[h p ci ts] = ttest(atanh(allRHO));
t = squeeze(ts.tstat); 

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end


%h(1:10) = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

mAR = mean(allRHO);
stdAR = std(allRHO); 
seAR = stdAR / sqrt(size(allRHO, 1))

xStart = -.2; dx = 0.05; times = xStart + (0:bins-1)*dx;
h(2, :) = times; 
figure()
%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'Set3')*0.75;
shadedErrorBar(times, mAR, seAR, {'Color',colors2use(1,:)}, 1); hold on; 

set(gca, 'ylim', [-.4 .2], 'xlim', [-.25 1.75]) % 

plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(times, hb, Linewidth=6); 

set(gca, 'Fontsize', 30)

exportgraphics(gcf, ['myP.png'], 'Resolution',150)



%% Correlate different trial-level metrics at ALL TIME POINTS NEW TIMES
clear, clc


minTrN = 8; 

printClust = 1; 
print1Clust = 0; 

paths = load_paths_EXT;

%c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
%c2 = 'trlSTA_HPC_CE_1-44_1_0_500-50'; 

c1 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 
c2 = 'trlCTX_TMP_CE_1-44_1_0_500-50'; 

%c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
%c2 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 

sub2exc = [37]; 
[cond1 cond2 ] = determine_conds_EXT(c1, c2, paths); 


nTimepoints = 51; 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins, bins); 
nTrials = zeros(50, 1); 

f2s = [c1(4:14) '_' c2(4:14)]; 

for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2)
        [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
        rsa2T1 = rsa2T1(i1, :); 
        rsa2T1IDs = rsa2T1IDs(i1,:); 

        rsa2T2 = rsa2T2(i2, :); 
        rsa2T2IDs = rsa2T2IDs(i2,:); 

        
        % % % % z-score both trial-level metrics separately for CS+ and CS-
        % rsa2TCSp = rsa2T1(rsa2T1IDs(:, 8) == 1,:); 
        % rsa2TCSm = rsa2T1(rsa2T1IDs(:, 8) == 0, :); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2T1(rsa2T1IDs(:, 8) == 1,:) = rsa2TCSp; 
        % rsa2T1(rsa2T1IDs(:, 8) == 0,:) = rsa2TCSm; 
        % 
        % 
        % rsa2TCSp = rsa2T2(rsa2T2IDs(:, 8) == 1,:); 
        % rsa2TCSm = rsa2T2(rsa2T2IDs(:, 8) == 0,:); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2T2(rsa2T2IDs(:, 8) == 1,:) = rsa2TCSp; 
        % rsa2T2(rsa2T2IDs(:, 8) == 0,:) = rsa2TCSm; 


        nTrials(subji, :) = size(rsa2T1, 1); 
    
        
         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 

            for timej = 1:bins
                timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
       
                allRHO(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
            end
        end
    end
end


sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc = unique(sub2exc3); 
allRHO(sub2exc, :, :) = []; 

nSubj = size(allRHO, 1); 
disp(['Number of subjects: ' num2str(nSubj)]); 

[h p ci ts] = ttest(atanh(allRHO), 0,'Alpha',0.05); 
%[h2 p ci ts] = ttest(atanh(allRHO), 0,'Alpha',0.01); 

h = squeeze(h); t = squeeze(ts.tstat); 
%h2 = squeeze(h2); 

h(1:5, :) = 0;
h(:, 1:5) = 0;




clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end

if ~printClust
    h = zeros(size(h)); 
end
if print1Clust
    h(clustinfo.PixelIdxList{id}) = 1; 
end

%h = zeros(size(h)); 
% h2(1:22, :) = 0;
% h2(:, 1:22) = 0;

 %h(clustinfo.PixelIdxList{2}) = 1; 
 % h(clustinfo.PixelIdxList{1}) = 1; 
% h(clustinfo.PixelIdxList{3}) = 1; 

d2p = squeeze(mean(allRHO)); 
figure(); colormap (brewermap([100], '*Spectral'))
contourf(myresizem(t, 10), 50, 'linecolor', 'none'); axis square; hold on; 
contour( myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
% contour( myresizem(h2, 10), 1, 'k:', 'LineWidth', 4);
%contour( myresizem(h, 10), 1, 'k:', 'Color', [0, 0, 0], 'LineWidth', 5); % % DASHED OUTLINE
plot([55 55],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(get(gca,'xlim'), [55 55],'k:', 'linewidth', 4); hold on; 
set(gca, 'xlim', [5 400], 'ylim', [5 400], 'FontSize', 24)
set(gca, 'xTick', [55 155 255 355], 'XTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'yTick', [55 155 255 355], 'YTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'clim', [-4 4])

exportgraphics(gcf, 'myP.png', 'Resolution', 300);



%% PERMUTATIONS
clearvars -except max_clust_obs cond1 cond2 sub2exc f2s
clc

% sub2exc inherited from previous code block

nPerm = 1000;
t4P = 6:40; 

nTimepoints = length(t4P); 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );


paths = load_paths_EXT; 


tic 

for permi = 1:nPerm

    clearvars -except nTimepoints win_width mf bins cond1 cond2 permi nPerm max_clust_perm max_clust_obs t4P sub2exc f2s paths

    progress_in_console(permi)

    allRHOP = zeros(50, bins, bins); 
    parfor subji = 1:50
    
        rsa2T1 = cond1{subji, 1}; 
        rsa2T1IDs = cond1{subji, 2}; 
    
        rsa2T2 = cond2{subji, 1}; 
        rsa2T2IDs = cond2{subji, 2}; 
        
        if ~isempty(rsa2T1) & ~isempty(rsa2T2)
            [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
            rsa2T1 = rsa2T1(i1, t4P); 
            rsa2T1IDs = rsa2T1IDs(i1,:); 
    
            rsa2T2 = rsa2T2(i2, t4P); 
            rsa2T2IDs = rsa2T2IDs(i2,:); 

            ids4perm = randperm(size(rsa2T2, 1)); 
            rsa2T2 = rsa2T2(ids4perm, :); 
            rsa2T2IDs = rsa2T2IDs(ids4perm, :); 
        
            
             for timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 
    
                for timej = 1:bins
                    timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                    rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
           
                    allRHOP(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
                end
            end
        end
    end
    
    allRHOP(sub2exc, :, :) = []; 
    
    [h p ci ts] = ttest(atanh(allRHOP));
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, : ) = allSTs(id); 
    else
        max_clust_perm(permi, : ) = 0; 
    end
    
    
end    

allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
%allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

f2s = [f2s '_' num2str(nPerm), 'p'];
save([paths.results.rsa_perm f2s], 'max_clust_perm', 'p', 'nPerm', 'max_clust_obs')


toc






%% Correlate different trial-level metrics at ALL TIME POINTS NEW TIMES SEPARATELY FOR EACH CONDITION
clear, clc

c2u1 = 1; %cs++ cs+- cs--
c2u2 = 1; % use the same in both c2u1 and c2u2 or different for OR

minTrN = 8; 

printClust = 1; 
print1Clust = 0; 

paths = load_paths_EXT;

c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_HPC_CE_1-44_1_0_500-50'; 

%c1 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 
%c2 = 'trlCTX_AMY_CE_1-44_1_0_500-50'; 

%c1 = 'trlSTA_TMP_CE_1-44_1_0_500-50'; 
%c2 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 

sub2exc = [37]; 
[cond1 cond2 ] = determine_conds_EXT(c1, c2, paths); 


nTimepoints = 51; 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins, bins); 
nTrials = zeros(50, 1); 

f2s = [c1(4:14) '_' c2(4:14)]; 

for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2)

        idh2uC1 = rsa2T1IDs(:, 6) == c2u1 | rsa2T1IDs(:, 6) == c2u2; 
        idh2uC2 = rsa2T2IDs(:, 6) == c2u1 | rsa2T2IDs(:, 6) == c2u2; 
        rsa2T1IDs = rsa2T1IDs(idh2uC1, :); 
        rsa2T2IDs = rsa2T2IDs(idh2uC2, :); 

        [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
        rsa2T1 = rsa2T1(i1, :); 
        rsa2T1IDs = rsa2T1IDs(i1,:); 

        rsa2T2 = rsa2T2(i2, :); 
        rsa2T2IDs = rsa2T2IDs(i2,:); 

       

        nTrials(subji, :) = size(rsa2T1, 1); 
    
        if ~isempty(rsa2T1) & ~isempty(rsa2T2) 
         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 

            for timej = 1:bins
                timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
       
                allRHO(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
            end
         end
        end
    end
end


sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc = unique(sub2exc3); 
allRHO(sub2exc, :, :) = []; 

nSubj = size(allRHO, 1); 
disp(['Number of subjects: ' num2str(nSubj)]); 

[h p ci ts] = ttest(atanh(allRHO));
h = squeeze(h); t = squeeze(ts.tstat); 

h(1:5, :) = 0;
h(:, 1:5) = 0;




clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end

if ~printClust
    h = zeros(size(h)); 
end
if print1Clust
    h(clustinfo.PixelIdxList{id}) = 1; 
end

%h = zeros(size(h)); 
% h(clustinfo.PixelIdxList{1}) = 1; 
% h(clustinfo.PixelIdxList{3}) = 1; 

d2p = squeeze(mean(allRHO)); 
figure(); colormap (brewermap([100], '*Spectral'))
contourf(myresizem(t, 10), 50, 'linecolor', 'none'); axis square; hold on; 
contour( myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
%contour( myresizem(h, 10), 1, 'k:', 'Color', [0, 0, 0], 'LineWidth', 5); % % DASHED OUTLINE
plot([55 55],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(get(gca,'xlim'), [55 55],'k:', 'linewidth', 4); hold on; 
set(gca, 'xlim', [5 400], 'ylim', [5 400], 'FontSize', 24)
set(gca, 'xTick', [55 155 255 355], 'XTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'yTick', [55 155 255 355], 'YTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'clim', [-4 4])

exportgraphics(gcf, 'myP.png', 'Resolution', 300);


%% PERMUTATIONS for separate in each condition analysis (first run previous cell)
clearvars -except max_clust_obs cond1 cond2 sub2exc f2s c2u1 c2u2
clc

% sub2exc inherited from previous code block

nPerm = 10;
t4P = 6:40; 

nTimepoints = length(t4P); 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );


paths = load_paths_EXT; 


tic 

for permi = 1:nPerm

    clearvars -except nTimepoints win_width mf bins cond1 cond2 permi nPerm max_clust_perm max_clust_obs t4P sub2exc f2s paths c2u1 c2u2

    progress_in_console(permi)

    allRHOP = zeros(50, bins, bins); 
    for subji = 1:50
    
        rsa2T1 = cond1{subji, 1}; 
        rsa2T1IDs = cond1{subji, 2}; 
    
        rsa2T2 = cond2{subji, 1}; 
        rsa2T2IDs = cond2{subji, 2}; 

        if ~isempty(rsa2T1) & ~isempty(rsa2T2)
        
        idh2uC1 = rsa2T1IDs(:, 6) == c2u1 | rsa2T1IDs(:, 6) == c2u2; 
        idh2uC2 = rsa2T2IDs(:, 6) == c2u1 | rsa2T2IDs(:, 6) == c2u2; 
        rsa2T1IDs = rsa2T1IDs(idh2uC1, :); 
        rsa2T2IDs = rsa2T2IDs(idh2uC2, :); 
        
        if ~isempty(rsa2T1IDs) & ~isempty(rsa2T2IDs)
            [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
            rsa2T1 = rsa2T1(i1, t4P); 
            rsa2T1IDs = rsa2T1IDs(i1,:); 
    
            rsa2T2 = rsa2T2(i2, t4P); 
            rsa2T2IDs = rsa2T2IDs(i2,:); 

            ids4perm = randperm(size(rsa2T2, 1)); 
            rsa2T2 = rsa2T2(ids4perm, :); 
            rsa2T2IDs = rsa2T2IDs(ids4perm, :); 
        
            if ~isempty(rsa2T1IDs) & ~isempty(rsa2T2IDs)
             parfor timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 
    
                for timej = 1:bins
                    timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                    rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
           
                    allRHOP(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
                end
             end
            end
        end
        end
    end
    
    allRHOP(sub2exc, :, :) = []; 
    
    [h p ci ts] = ttest(atanh(allRHOP));
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, : ) = allSTs(id); 
    else
        max_clust_perm(permi, : ) = 0; 
    end
    
    
end    

allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
%allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

f2s = [f2s '_' num2str(c2u1) '_' num2str(c2u2) '_' num2str(nPerm), 'p']
save([paths.results.rsa_perm f2s], 'max_clust_perm', 'p', 'nPerm', 'max_clust_obs')


toc



%% Correlate different trial-level metrics at ALL TIME POINTS NEW TIMES SEPARATELY CONTRAST BETWEEN CONDITIONS
clear, clc

c2u1 = 1; %cs++ cs+- cs--
c2u2 = 2; 
c2u3 = 2; % % the structure is c2u1 vs c2u2 OR C2u3

minTrN = 8; 

printClust = 1; 
print1Clust = 0; 

paths = load_paths_EXT;

%c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
%c2 = 'trlSTA_HPC_CE_1-44_1_0_500-50'; 

c1 = 'trlCTX_TMP_CE_1-44_1_0_500-50'; 
c2 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 

%c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
%c2 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 

sub2exc = [37]; 
[cond1 cond2 ] = determine_conds_EXT(c1, c2, paths); 


nTimepoints = 51; 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO1 = zeros(50, bins, bins); 
allRHO2 = zeros(50, bins, bins); 
nTrials = zeros(50, 1); 

f2s = [c1(4:14) '_' c2(4:14)]; 

for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2)

        idh2uC1 = rsa2T1IDs(:, 6) == c2u1 ; 
        idh2uC2 = rsa2T2IDs(:, 6) == c2u1 ; 
        rsa2T1IDs = rsa2T1IDs(idh2uC1, :); 
        rsa2T2IDs = rsa2T2IDs(idh2uC2, :); 

        [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
        rsa2T1 = rsa2T1(i1, :); 
        rsa2T1IDs = rsa2T1IDs(i1,:); 

        rsa2T2 = rsa2T2(i2, :); 
        rsa2T2IDs = rsa2T2IDs(i2,:); 

       

        nTrials(subji, :) = size(rsa2T1, 1); 
    
        if ~isempty(rsa2T1) & ~isempty(rsa2T2) 
         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 

            for timej = 1:bins
                timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
       
                allRHO1(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
            end
         end
        end
    end
end


sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc1 = unique(sub2exc3); 




for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2)

        idh2uC1 = rsa2T1IDs(:, 6) == c2u2 | rsa2T1IDs(:, 6) == c2u3; 
        idh2uC2 = rsa2T2IDs(:, 6) == c2u2 | rsa2T2IDs(:, 6) == c2u3; 
        rsa2T1IDs = rsa2T1IDs(idh2uC1, :); 
        rsa2T2IDs = rsa2T2IDs(idh2uC2, :); 

        [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
        rsa2T1 = rsa2T1(i1, :); 
        rsa2T1IDs = rsa2T1IDs(i1,:); 

        rsa2T2 = rsa2T2(i2, :); 
        rsa2T2IDs = rsa2T2IDs(i2,:); 

        nTrials(subji, :) = size(rsa2T1, 1); 
    
        if ~isempty(rsa2T1) & ~isempty(rsa2T2) 
         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 

            for timej = 1:bins
                timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
       
                allRHO2(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
            end
         end
        end
    end
end


sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc2 = unique(sub2exc3); 

sub2exc4 = union(sub2exc1, sub2exc2); 

allRHO1(sub2exc4, :, :) = []; 
allRHO2(sub2exc4, :, :) = []; 

nSubj = size(allRHO1, 1); 
disp(['Number of subjects: ' num2str(nSubj)]); 

[h p ci ts] = ttest(atanh(allRHO1), atanh(allRHO2));
h = squeeze(h); t = squeeze(ts.tstat); 

h(1:5, :) = 0;
h(:, 1:5) = 0;



clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end

if ~printClust
    h = zeros(size(h)); 
end
if print1Clust
    h(clustinfo.PixelIdxList{id}) = 1; 
end

% h = zeros(size(h)); 
% h(clustinfo.PixelIdxList{1}) = 1; 
% h(clustinfo.PixelIdxList{3}) = 1; 

%d2p = squeeze(mean(allRHO)); 
figure(); colormap (brewermap([100], '*Spectral'))
contourf(myresizem(t, 10), 50, 'linecolor', 'none'); axis square; hold on; 
contour( myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
%contour( myresizem(h, 10), 1, 'k:', 'Color', [0, 0, 0], 'LineWidth', 5); % % DASHED OUTLINE
plot([55 55],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(get(gca,'xlim'), [55 55],'k:', 'linewidth', 4); hold on; 
set(gca, 'xlim', [5 400], 'ylim', [5 400], 'FontSize', 24)
set(gca, 'xTick', [55 155 255 355], 'XTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'yTick', [55 155 255 355], 'YTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'clim', [-4 4])

exportgraphics(gcf, 'myP.png', 'Resolution', 300);

%% permutations (shuffling condition labels in every subject)

clc

% sub2exc inherited from previous code block
nPerm = 1000;
tic
for permi = 1:nPerm

    progress_in_console(permi)
    clear allRHO1p allRHO2p
    for subji = 1:size(allRHO1, 1)
        if rand>.5
           allRHO1p(subji, :, :) =  allRHO2(subji, :, :) ; 
           allRHO2p(subji, :, :) =  allRHO1(subji, :, :) ; 
        else
           allRHO1p(subji, :, :) =  allRHO1(subji, :, :) ; 
           allRHO2p(subji, :, :) =  allRHO2(subji, :, :) ; 
        end
    end

    allRHO1p = allRHO1p(:, 6:42, 6:42); 
    allRHO2p = allRHO2p(:, 6:42, 6:42); 

    [h p ci ts] = ttest(atanh(allRHO1p), atanh(allRHO2p));
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, : ) = allSTs(id); 
    else
        max_clust_perm(permi, : ) = 0; 
    end
    
    
end    

allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
%allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

f2s = [f2s '_' num2str(c2u1) '_' num2str(c2u2) '_' num2str(c2u3) '_' num2str(nPerm), 'p']
save([paths.results.rsa_perm f2s], 'max_clust_perm', 'p', 'nPerm', 'max_clust_obs')


toc



%% correlate CTX SPE in PFC across SUBJECTS with REINST in AMY AND TMP (across subjects)


clear, clc

trltype = 3;
minTrN = 8; 
tP = 21:28;  %21:28 PFC effect %19:31 PFC effec ACQvsEXT

remOutliers = 0; 
sub2exc = [37]'; 


c1 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_TMP_CAT_1-44_1_0_500-50'; 
c3 = 'trlSTA_TMP_CET_1-44_1_0_500-50'; 


paths = load_paths_EXT;


[cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
[cond1 cond3] = determine_conds_EXT(c1, c3, paths); 

rsa2TT1 = zeros(50, 1); 
rsa2TT2 = zeros(50, 1); 
rsa2TT3 = zeros(50, 1); 


for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 

    rsa2T3 = cond3{subji, 1}; 
    rsa2T3IDs = cond3{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2) & ~isempty(rsa2T3)
        
        if trltype ~= 0
            rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == trltype, tP); 
            rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, tP); 
            rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == trltype, tP); 
        else 
            rsa2T1 = rsa2T1(:, tP); 
            rsa2T2 = rsa2T2(:, tP); 
            rsa2T3 = rsa2T3(:, tP); 
        end
        

        nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 

        rsa2TT1(subji, :) = mean(rsa2T1, 'all', 'omitnan'); 
        rsa2TT2(subji, :) = mean(rsa2T2, 'all', 'omitnan'); 
        rsa2TT3(subji, :) = mean(rsa2T3, 'all', 'omitnan'); 

    end
end

sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc; sub2exc2]; 
sub2exc = unique(sub2exc3); 

rsa2TT1(sub2exc, :, :) = []; 
rsa2TT2(sub2exc, :, :) = []; 
rsa2TT3(sub2exc, :, :) = []; 

rsa2TT1(rsa2TT1==0) = []; 
rsa2TT2(rsa2TT2==0) = []; 
rsa2TT3(rsa2TT3==0) = []; 

rsa2TT4 = rsa2TT2 - rsa2TT3;  %ACQ minus EXT
%rsa2TT4 = rsa2TT3; 

if remOutliers
    %[B, tF1] = rmoutliers(rsa2TT1, 'percentile', [5 95]); 
    %[B, tF2] = rmoutliers(rsa2TT4, 'percentile', [5 95]); 
    [B, tF1] = rmoutliers(rsa2TT1, 'mean'); 
    [B, tF2] = rmoutliers(rsa2TT4, 'mean'); 
    tF3 = tF1 | tF2; 
    nOut(subji, :) = sum(tF3); 
    rsa2TT1(tF3) = []; 
    rsa2TT4(tF3) = []; 
end

 
[rho p] = corr(rsa2TT1, rsa2TT4, 'type', 's');

nSubj = length(rsa2TT1); 
disp(['Number of subjects: ' num2str(nSubj)]); 
scatter(rsa2TT1, rsa2TT4, 1400, '.'); hold on; 
scatter(rsa2TT1, rsa2TT4, 400, 'ko'); hold on; 
pFit = polyfit(rsa2TT1, rsa2TT4, 1);
m = pFit(1); % slope
b = pFit(2); % intercept
line([min(rsa2TT1)-.03 max(rsa2TT1)+.03], [m*min(rsa2TT1)+b m*max(rsa2TT1)+b], 'color', 'r', linewidth=3);
%title (['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]);

set(gca, Fontsize=28)
set(gca, xlim=[-.035 .05], ylim=[-.06 .06])
exportgraphics(gcf, 'myP.png', 'Resolution', 300);
disp(['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]); 



%% correlate ITEM STABILITY WITH with REINST in AMY AND TMP (across subjects)


clear, clc

trltype = 3;
minTrN = 8; 
%tP = 21:28;  %21:28 PFC effect %19:31 PFC effec ACQvsEXT
tP = 29:33; %29:33; %18:25;  %21:28 


remOutliers = 0; 
sub2exc = [37]'; 


c1 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_TMP_CAT_1-44_1_0_500-50'; 
c3 = 'trlSTA_TMP_CET_1-44_1_0_500-50'; 


paths = load_paths_EXT;


[cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
[cond1 cond3] = determine_conds_EXT(c1, c3, paths); 

rsa2TT1 = zeros(50, 1); 
rsa2TT2 = zeros(50, 1); 
rsa2TT3 = zeros(50, 1); 


for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 

    rsa2T3 = cond3{subji, 1}; 
    rsa2T3IDs = cond3{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2) & ~isempty(rsa2T3)
        
        if trltype ~= 0
            rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == trltype, tP); 
            rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, tP); 
            rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == trltype, tP); 
            % rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == 1 | rsa2T1IDs(:, 6) == 3, tP); 
            % rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == 1 | rsa2T2IDs(:, 6) == 3, tP); 
            % rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == 1 | rsa2T3IDs(:, 6) == 3, tP); 
        else 
            rsa2T1 = rsa2T1(:, tP); 
            rsa2T2 = rsa2T2(:, tP); 
            rsa2T3 = rsa2T3(:, tP); 
        end
        

        nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 

        rsa2TT1(subji, :) = mean(rsa2T1, 'all', 'omitnan'); 
        rsa2TT2(subji, :) = mean(rsa2T2, 'all', 'omitnan'); 
        rsa2TT3(subji, :) = mean(rsa2T3, 'all', 'omitnan'); 

    end
end

sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc; sub2exc2]; 
sub2exc = unique(sub2exc3); 

rsa2TT1(sub2exc, :, :) = []; 
rsa2TT2(sub2exc, :, :) = []; 
rsa2TT3(sub2exc, :, :) = []; 

rsa2TT1(rsa2TT1==0) = []; 
rsa2TT2(rsa2TT2==0) = []; 
rsa2TT3(rsa2TT3==0) = []; 

rsa2TT4 = rsa2TT2 - rsa2TT3;  %ACQ minus EXT
%rsa2TT4 = rsa2TT3; 

if remOutliers
    %[B, tF1] = rmoutliers(rsa2TT1, 'percentile', [5 95]); 
    %[B, tF2] = rmoutliers(rsa2TT4, 'percentile', [5 95]); 
    [B, tF1] = rmoutliers(rsa2TT1, 'mean'); 
    [B, tF2] = rmoutliers(rsa2TT4, 'mean'); 
    tF3 = tF1 | tF2; 
    nOut(subji, :) = sum(tF3); 
    rsa2TT1(tF3) = []; 
    rsa2TT4(tF3) = []; 
end

 
[rho p] = corr(rsa2TT1, rsa2TT4, 'type', 's');

nSubj = length(rsa2TT1); 
disp(['Number of subjects: ' num2str(nSubj)]); 
scatter(rsa2TT1, rsa2TT4, 1400, '.'); hold on; 
scatter(rsa2TT1, rsa2TT4, 400, 'ko'); hold on; 
pFit = polyfit(rsa2TT1, rsa2TT4, 1);
m = pFit(1); % slope
b = pFit(2); % intercept
line([min(rsa2TT1)-.03 max(rsa2TT1)+.03], [m*min(rsa2TT1)+b m*max(rsa2TT1)+b], 'color', 'r', linewidth=3);
%title (['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]);

set(gca, Fontsize=28)
set(gca, xlim=[-.035 .05], ylim=[-.06 .06])
exportgraphics(gcf, 'myP.png', 'Resolution', 300);
disp(['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]); 


%% Correlate amygdala POWER IN CLUSTER with different trial-level metrics (ALL TIME POINTS)
clear, clc

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_1-44_1_0_500-50'])



load ([paths.results.trial_based 'trlSTA_AMY_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_1-44_1_0_500-50'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_1-44_1_0_500-50'])

sub2exc = [37]; 
minTr = 5; 



if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end



win_width = 10; 
mf = 1; 
nTimepoints = 51; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins); 

for subji = 1:50

    amyPOW1 = allPOWAMY{subji, 1}; 
    amyPOWIDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 
    
    if ~isempty(amyPOW1) & ~isempty(rsa2T)
        [i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW1(i2, :); 
        amyPOWIDsh = amyPOWIDs(i2,:); 
    
        nanIds = isnan(amyPOW); 
        amyPOW(nanIds,:) = []; 
        rsa2T(nanIds,:) = []; 
        amyPOWIDsh(nanIds,:) = []; 
        rsa2TIDs(nanIds,:) = []; 


        % % % z-score amygdala power separately for CS+ and CS-
        amyPOWCSp = amyPOW(amyPOWIDsh(:, 8) == 1); 
        amyPOWCSm = amyPOW(amyPOWIDsh(:, 8) == 0); 
        amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
        amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
        amyPOW(amyPOWIDsh(:, 8) == 1) = amyPOWCSp; 
        amyPOW(amyPOWIDsh(:, 8) == 0) = amyPOWCSm; 
        % % % z-score RSA2T metric separately for CS+ and CS-
        % rsa2TCSp = rsa2TT(rsa2TIDs(:, 8) == 1); 
        % rsa2TCSm = rsa2TT(rsa2TIDs(:, 8) == 0); 
        % rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        % rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        % rsa2TT(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        % rsa2TT(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nTrials (subji, : ) = length(amyPOW); 
        



         for timei = 1:bins 
            %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
            rsa2TT = mean(rsa2T(:, timeBins), 2);     
    
            % [B, tF1] = rmoutliers(rsa2TT, 'percentiles', [5 95]); 
            % [B, tF2] = rmoutliers(amyPOW, 'percentiles', [5 95]); 
            % [B, tF1] = rmoutliers(rsa2TT, 'mean'); 
            % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
            % tF3 = tF1 | tF2; 
            % nOut(subji, timei) = sum(tF3); 
            % rsa2TT(tF3,:) = []; 
            % amyPOWH = amyPOW; 
            % amyPOWH(tF3,:) = []; 

            if length(amyPOW) > minTr % at least five trials for the correlation analysis
                allRHO(subji, timei, :) = corr(amyPOW, rsa2TT, 'type', 's');
            
            else
            
                allRHO(subji, timei, :) = nan; 
            end
    
         end

    end
end

allRHO(sub2exc, :) = []; 

allRHO = allRHO(any(allRHO,2),:);
[h p ci ts] = ttest(atanh(allRHO));
t = squeeze(ts.tstat); 

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end


%h(1:10) = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

mAR = mean(allRHO);
stdAR = std(allRHO); 
seAR = stdAR / sqrt(size(allRHO, 1))

xStart = -.2; dx = 0.05; times = xStart + (0:bins-1)*dx;
h(2, :) = times; 
figure()
%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'Set3')*0.75;
shadedErrorBar(times, mAR, seAR, {'Color',colors2use(1,:)}, 1); hold on; 

set(gca, 'ylim', [-.4 .2], 'xlim', [-.25 1.75]) % 

plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(times, hb, Linewidth=6); 

set(gca, 'Fontsize', 30)

exportgraphics(gcf, ['myP.png'], 'Resolution',150)






%% correlate THETA POWER WITH with REINST in AMY AND TMP (across subjects)


clear, clc

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])

trltype = 3;
minTrN = 8; 
tP = 29:33;  %21:28 PFC effect %19:31 PFC effec ACQvsEXT
%tP = 18:25; %29:33; %18:25;  %21:28 


remOutliers = 0; 
sub2exc = [37]'; 

c2 = 'trlSTA_AMY_CAT_1-44_1_0_500-50'; 
c3 = 'trlSTA_AMY_CET_1-44_1_0_500-50'; 


paths = load_paths_EXT;


[cond2 cond3] = determine_conds_EXT(c2, c3, paths); 
%[cond1 cond3] = determine_conds_EXT(c1, c3, paths); 

rsa2TT1 = zeros(50, 1); 
rsa2TT2 = zeros(50, 1); 
rsa2TT3 = zeros(50, 1); 


for subji = 1:50

    rsa2T1 = allPOWAMY{subji, 1}; 
    rsa2T1IDs = double(string(allPOWAMY{subji, 2})); 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 

    rsa2T3 = cond3{subji, 1}; 
    rsa2T3IDs = cond3{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2) & ~isempty(rsa2T3)
        
        if trltype ~= 0
            %rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == trltype, :); 
            %rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, tP); 
            %rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == trltype, tP); 
             rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == 2 | rsa2T1IDs(:, 6) == 3, :); 
             rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == 2 | rsa2T2IDs(:, 6) == 3, tP); 
             rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == 2 | rsa2T3IDs(:, 6) == 3, tP); 
        else 
            rsa2T1 = rsa2T1; 
            rsa2T2 = rsa2T2(:, tP); 
            rsa2T3 = rsa2T3(:, tP); 
        end
        

        nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 

        rsa2TT1(subji, :) = mean(rsa2T1, 'all', 'omitnan'); 
        rsa2TT2(subji, :) = mean(rsa2T2, 'all', 'omitnan'); 
        rsa2TT3(subji, :) = mean(rsa2T3, 'all', 'omitnan'); 

    end
end

sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc; sub2exc2]; 
sub2exc = unique(sub2exc3); 

rsa2TT1(sub2exc, :, :) = []; 
rsa2TT2(sub2exc, :, :) = []; 
rsa2TT3(sub2exc, :, :) = []; 

rsa2TT1(rsa2TT1==0) = []; 
rsa2TT2(rsa2TT2==0) = []; 
rsa2TT3(rsa2TT3==0) = []; 

rsa2TT4 = rsa2TT2 - rsa2TT3;  %ACQ minus EXT
%rsa2TT4 = rsa2TT3; 

if remOutliers
    %[B, tF1] = rmoutliers(rsa2TT1, 'percentile', [5 95]); 
    %[B, tF2] = rmoutliers(rsa2TT4, 'percentile', [5 95]); 
    [B, tF1] = rmoutliers(rsa2TT1, 'mean'); 
    [B, tF2] = rmoutliers(rsa2TT4, 'mean'); 
    tF3 = tF1 | tF2; 
    nOut(subji, :) = sum(tF3); 
    rsa2TT1(tF3) = []; 
    rsa2TT4(tF3) = []; 
end

 
[rho p] = corr(rsa2TT1, rsa2TT4, 'type', 's');

nSubj = length(rsa2TT1); 
disp(['Number of subjects: ' num2str(nSubj)]); 
scatter(rsa2TT1, rsa2TT4, 1400, '.'); hold on; 
scatter(rsa2TT1, rsa2TT4, 400, 'ko'); hold on; 
pFit = polyfit(rsa2TT1, rsa2TT4, 1);
m = pFit(1); % slope
b = pFit(2); % intercept
line([min(rsa2TT1)-.03 max(rsa2TT1)+.03], [m*min(rsa2TT1)+b m*max(rsa2TT1)+b], 'color', 'r', linewidth=3);
%title (['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]);

set(gca, Fontsize=28)
set(gca, xlim=[-.2 .2], ylim=[-.2 .2])
exportgraphics(gcf, 'myP.png', 'Resolution', 300);
disp(['Rho: ' num2str(rho, 3) ' p: ' num2str(p, 3)]); 

%% correlate across SUBJECTS and not trials : ONE time period PFC, ALL TIME PERIODS IN TMP

clear, clc


trltype = 2;
tP = 21:28; % PFC CTX cluster; 
minTrN = 8; 
win_width = 8; 


paths = load_paths_EXT;

sub2exc =[37]'; 

c1 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_TMP_CAT_1-44_1_0_500-50'; 
c3 = 'trlSTA_TMP_CET_1-44_1_0_500-50'; 

[cond1 cond2 ] = determine_conds_EXT(c1, c2, paths); 
[cond1 cond3 ] = determine_conds_EXT(c1, c3, paths); 

rsa2TT1 = zeros(50, 1); 
rsa2TT2 = zeros(50, 1); 
rsa2TT3 = zeros(50, 1); 

nTimepoints = 51; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

rsa2TT1 = zeros(50, 1); 
rsa2TT2 = zeros(50, bins); 
rsa2TT3 = zeros(50, bins); 

for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 

    rsa2T3 = cond3{subji, 1}; 
    rsa2T3IDs = cond3{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2) & ~isempty(rsa2T3)
        
        if trltype ~= 0
            rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == trltype, :); 
            rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, :); 
            rsa2T3 = rsa2T3(rsa2T3IDs(:, 6) == trltype, :); 
        else 
            rsa2T1 = rsa2T1; 
            rsa2T2 = rsa2T2; 
            rsa2T3 = rsa2T3; 
        end
        
        nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 

        rsa2TT1(subji, :) = mean(mean(rsa2T1(:, tP), 2)); 

        for timej = 1:bins
            timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
            
            rsa2TT2(subji, timej) = mean(mean(rsa2T2(:, timeBinsj), 2)); 
            rsa2TT3(subji, timej) = mean(mean(rsa2T3(:, timeBinsj), 2)); 
   
        end
        
       

    end
end

sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc; sub2exc2]; 
sub2exc = unique(sub2exc3); 

rsa2TT1(sub2exc, :, :) = []; 
rsa2TT2(sub2exc, :, :) = []; 
rsa2TT3(sub2exc, :, :) = []; 
 
rsa2TT1 = rsa2TT1(any(rsa2TT1,2),:);
rsa2TT2 = rsa2TT2(any(rsa2TT2,2),:);
rsa2TT3 = rsa2TT3(any(rsa2TT3,2),:);

rsa2TT4 = rsa2TT2 - rsa2TT3; 

for timei = 1:bins
    [rho(timei, :) p(timei, :)] = corr(rsa2TT1, rsa2TT4(:, timei), 'type', 's');
end
h = p<0.05; hb = double(h); hb(hb==0) = nan; hb(hb==1) = -0.8; 



clustinfo = bwconncomp(h);
clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(rho(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end



nSubj = length(rsa2TT1); 

xStart = -.2; dx = 0.05; times = xStart + (0:bins-1)*dx;

plot(times, rho, LineWidth=5); hold on; 
%plot(times, hb,'k:',  LineWidth=7)
plot(times, hb,':',  LineWidth=7)
set(gca, fontsize=28, xlim=[-.2 1.75], ylim=[-1 1])
plot(get(gca, 'xlim'), [0 0], 'k:', 'linewidth', 3);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 3);
exportgraphics(gcf, 'myP.png', 'Resolution', 300);

%% PEMUTATIONS 

clearvars -except max_clust_obs rsa2TT1 rsa2TT4 win_width nTimepoints mf bins

nPerm = 1000; 


for permi = 1:nPerm

    progress_in_console(permi)


    ids4perm = randperm(size(rsa2TT1, 1));
    rsa2TT1P = rsa2TT1(ids4perm, :); 
    
    for timei = 1:bins
        [rho(timei, :) p(timei, :)] = corr(rsa2TT1P, rsa2TT4(:, timei), 'type', 's');
    end
    
    
    h = p<0.05; 
    h(1:5) = 0; 
    clustinfo = bwconncomp(h);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(rho(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, :) = allSTs(id); 
    else
        max_clust_perm(permi, :) = 0; 
    end




end

allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
%allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm







%% COMPUTE REINSTATEMENT SEPARATELY FOR CS++ CS+- and CS-- 


clear, clc

roi = 'AMY'; 
aet = 1; % 0 = ACQ-EXT; 1 = ACQ; 2 = EXT; 
minTrN = 8; 
%tP = 21:28;  %21:28 PFC effect %19:31 PFC effec ACQvsEXT
%tP = 18:25; %TMP ITEM STABILITY
tP = 29:33; %AMY ITEM STABILITY 


allSUB2EXC =  []; 
for typei = 1:3

    clearvars -except allRSATT typei allSUB2EXC aet roi tP minTrN 
    
    
    trltype = typei;

    
    remOutliers = 0; 
    sub2exc = [37]'; 
    
        
    c1 = ['trlSTA_' roi '_CAT_1-44_1_0_500-50']; 
    c2 = ['trlSTA_' roi '_CET_1-44_1_0_500-50']; 

    
    
    paths = load_paths_EXT;
    
    
    [cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
    
    rsa2TT1 = zeros(50, 1); 
    rsa2TT2 = zeros(50, 1); 
    

    for subji = 1:50
    
        rsa2T1 = cond1{subji, 1}; 
        rsa2T1IDs = cond1{subji, 2}; 
    
        rsa2T2 = cond2{subji, 1}; 
        rsa2T2IDs = cond2{subji, 2}; 
    
        
        if ~isempty(rsa2T1) & ~isempty(rsa2T2) 
            
            if trltype ~= 0
                rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) == trltype, tP); 
                rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, tP); 
                
            else 
                rsa2T1 = rsa2T1(:, tP); 
                rsa2T2 = rsa2T2(:, tP); 
            end
            
    
            %nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 
            nTrials(subji, :) = min([size(rsa2T1, 1),  size(rsa2T2, 1)]); 
    
            rsa2TT1(subji, :) = mean(rsa2T1, 'all', 'omitnan'); 
            rsa2TT2(subji, :) = mean(rsa2T2, 'all', 'omitnan'); 
            
    
        end
    end
    
    sub2exc2 = find(nTrials < minTrN); 
    sub2exc3 = [sub2exc; sub2exc2]; 
    sub2exc = unique(sub2exc3); 

    allSUB2EXC = [allSUB2EXC sub2exc']; 
    
    % rsa2TT1(sub2exc, :, :) = []; 
    % rsa2TT2(sub2exc, :, :) = []; 
    % rsa2TT3(sub2exc, :, :) = []; 
    % 
    % rsa2TT1(rsa2TT1==0) = []; 
    % rsa2TT2(rsa2TT2==0) = []; 
    % rsa2TT3(rsa2TT3==0) = []; 
    
    %rsa2TT3 = rsa2TT1 - rsa2TT2;  %ACQ minus EXT
    if aet == 1
        rsa2TT3 = rsa2TT1; 
    elseif aet == 2
        rsa2TT3 = rsa2TT2; 
    elseif aet == 0
        rsa2TT3 = rsa2TT1 - rsa2TT2;  %ACQ minus EXT
    end

    if remOutliers
        %[B, tF1] = rmoutliers(rsa2TT1, 'percentile', [5 95]); 
        %[B, tF2] = rmoutliers(rsa2TT4, 'percentile', [5 95]); 
        [B, tF1] = rmoutliers(rsa2TT1, 'mean'); 
        [B, tF2] = rmoutliers(rsa2TT4, 'mean'); 
        tF3 = tF1 | tF2; 
        nOut(subji, :) = sum(tF3); 
        rsa2TT1(tF3) = []; 
        rsa2TT4(tF3) = []; 
    end
    
    
    allRSATT{typei, :} = rsa2TT3; 
end

allSUB2EXC = unique(allSUB2EXC); 

for typei = 1:3

    allRSAH = allRSATT{typei, :} ; 
    allRSAH(allSUB2EXC, :, :) = []; 
    allRSAH(allRSAH==0) = [];
    allRSATT{typei, :}  = allRSAH; 
end

%%ANOVA REPEATED MEASURES

clc 
data = [allRSATT{1}, allRSATT{2}, allRSATT{3}];
nSubj = size(allRSATT{1}, 1); 
d4anova = data(:);
d4anova = d4anova(:); 
d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];

[p f] = RMAOV1(d4anova);

boxplot(data)

%%

c2u = 2; 
c1 = data_1(:, c2u); 
c2 = data_2(:, c2u); 

[h p ci ts] = ttest(c1, c2); 
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ';' ' p = ' num2str(p, 3)]);
    


%% plot two bar

ylim = [-.075 .15];
xlim = [0 4];

%data.data = [plvCSPPE plvCSMME]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data, 1);
std_S = std(data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1 2 3], data, 'k'); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',20);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2 3],'XTickLabel',{'   '},     'FontSize', 20, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
exportgraphics(gcf, 'myP.png', 'Resolution',150)

%%

[h p ci ts] = ttest(data(:, 1), data(:, 3))

%% Repeated-Measures ANOVA (Same as above)

% Combine your data into a matrix where each column corresponds to a condition
data = [allRSATT{1}, allRSATT{2}, allRSATT{3}];

% Create a table with valid variable names for each condition.
% (We use names like 'CSpp', 'CSpm', and 'CSmm' instead of 'CS++', etc.)
T = array2table(data, 'VariableNames', {'CSpp', 'CSpm', 'CSmm'});

% Create a within-subjects design table.
% The table must have as many rows as the number of repeated measures (here, 3 conditions).
withinDesign = table((1:3)', 'VariableNames', {'Condition'});
withinDesign.Condition = categorical(withinDesign.Condition);

% Fit the repeated measures model.
% The model specification 'CSpp-CSmm ~ 1' tells MATLAB to use the table columns from CSpp to CSmm.
rm = fitrm(T, 'CSpp-CSmm ~ 1', 'WithinDesign', withinDesign);

% Run the repeated-measures ANOVA.
% Specify the within-subject factor using the name you gave in withinDesign (here, 'Condition').
ranovatbl = ranova(rm, 'WithinModel', 'Condition');
disp('Repeated Measures ANOVA Table:');
disp(ranovatbl);

% (Optional) Do post-hoc pairwise comparisons with Bonferroni correction:
comp = multcompare(rm, 'Condition', 'ComparisonType', 'bonferroni');
disp('Pairwise Comparisons (Bonferroni corrected):');
disp(comp);


 
%% COMPUTE REINSTATEMENT SEPARATELY FOR CS++ CS+- and CS-- FOR ALL TIME POINTS


clear, clc

aet = 0; % 0 = ACQ-EXT; 1 = ACQ; 2 = EXT; 
roi = 'AMY';
minTrN = 8; 
remOutliers = 0; 
sub2exc = [37]'; 


allSUB2EXC =  []; 
nTimepoints = 51; 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );
allData = [];              
for timei = 1:bins 
    %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
    for typei = 1:3
    
        clearvars -except allRSATT typei allSUB2EXC timei allPs allFs timeBinsi mf win_width nTimepoints aet roi ...
            minTrN remOutliers sub2exc allData
        
        trltype = typei;
        
        tP = timeBinsi; %21:28;  %21:28 PFC effect %19:31 PFC effec ACQvsEXT
        %tP = 15:28;
        

        
        c1 = ['trlSTA_' roi '_CAT_1-44_1_0_500-50']; 
        c2 = ['trlSTA_' roi '_CET_1-44_1_0_500-50']; 
        
        
        paths = load_paths_EXT;
        
        
        [cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
        
        rsa2TT1 = zeros(50, 1); 
        rsa2TT2 = zeros(50, 1); 
        
    
        for subji = 1:50
        
            rsa2T1 = cond1{subji, 1}; 
            rsa2T1IDs = cond1{subji, 2}; 
        
            rsa2T2 = cond2{subji, 1}; 
            rsa2T2IDs = cond2{subji, 2}; 
        
            
            if ~isempty(rsa2T1) & ~isempty(rsa2T2) 
                
                if trltype ~= 0
                    rsa2T1= rsa2T1(rsa2T1IDs(:, 6) == trltype, tP); 
                    rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) == trltype, tP); 
                    
                else 
                    rsa2T1 = rsa2T1(:, tP); 
                    rsa2T2 = rsa2T2(:, tP); 
                end
                
        
                %nTrials(subji, :) = min([size(rsa2T1, 1), size(rsa2T2, 1),  size(rsa2T3, 1)]); 
                nTrials(subji, :) = min([size(rsa2T1, 1),  size(rsa2T2, 1)]); 
        
                rsa2TT1(subji, :) = mean(rsa2T1, 'all', 'omitnan'); 
                rsa2TT2(subji, :) = mean(rsa2T2, 'all', 'omitnan'); 
                
        
            end
        end
        
        sub2exc2 = find(nTrials < minTrN); 
        sub2exc3 = [sub2exc; sub2exc2]; 
        sub2exc = unique(sub2exc3); 
    
        allSUB2EXC = [allSUB2EXC sub2exc']; 
        
        % rsa2TT1(sub2exc, :, :) = []; 
        % rsa2TT2(sub2exc, :, :) = []; 
        % rsa2TT3(sub2exc, :, :) = []; 
        % 
        % rsa2TT1(rsa2TT1==0) = []; 
        % rsa2TT2(rsa2TT2==0) = []; 
        % rsa2TT3(rsa2TT3==0) = []; 
        
        %rsa2TT3 = rsa2TT1 - rsa2TT2;  %ACQ minus EXT
        if aet == 1
            rsa2TT3 = rsa2TT1; 
        elseif aet == 2
            rsa2TT3 = rsa2TT2; 
        elseif aet == 0
            rsa2TT3 = rsa2TT1 - rsa2TT2;  %ACQ minus EXT
        end
    
        if remOutliers
            %[B, tF1] = rmoutliers(rsa2TT1, 'percentile', [5 95]); 
            %[B, tF2] = rmoutliers(rsa2TT4, 'percentile', [5 95]); 
            [B, tF1] = rmoutliers(rsa2TT1, 'mean'); 
            [B, tF2] = rmoutliers(rsa2TT4, 'mean'); 
            tF3 = tF1 | tF2; 
            nOut(subji, :) = sum(tF3); 
            rsa2TT1(tF3) = []; 
            rsa2TT4(tF3) = []; 
        end
        
        
        allRSATT{typei, :} = rsa2TT3; 
        
    end
    
    allSUB2EXC = unique(allSUB2EXC); 
    
    for typei = 1:3
    
        allRSAH = allRSATT{typei, :} ; 
        allRSAH(allSUB2EXC, :, :) = []; 
        allRSAH(allRSAH==0) = [];
        allRSATT{typei, :}  = allRSAH; 
    end
    
    %%ANOVA REPEATED MEASURES
    
    clc 
    data = [allRSATT{1}, allRSATT{2}, allRSATT{3}];
    allData(:, :, timei) = data; 
    nSubj = size(allRSATT{1}, 1); 
    d4anova = data(:);
    d4anova = d4anova(:); 
    d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
    d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];
    
    [allPs(timei, :) allFs(timei, :)] = RMAOV1(d4anova);
    
    
end

figure()
times = -.2:.05:1.85;
plot (times, allFs); hold on; 
plot(times, allPs)
set(gca, 'xlim', [-.25 1.7]); 

%% plot tree conditions
clc 

mD = squeeze(mean(allData, 1)); 
stD = squeeze(std(allData, [], 1)); 
se = stD / sqrt(size(allData, 1)); 

fig_stuff=subplot(1,1,1);
cmap_default=fig_stuff.ColorOrder;
green= colormap(brewermap([1],'Greens'))
green = green*.9;
red = cmap_default(2,:);
yellow = cmap_default(3,:); 

times = -.2:.05:1.85;
boundedline(times, mD(1,:), se(1, :),'LineWidth', 2, 'cmap',red,'transparency',0.2,'alpha'); hold on; 
boundedline(times, mD(2,:), se(2, :),'LineWidth', 2, 'cmap',yellow,'transparency',0.2,'alpha');
boundedline(times, mD(3,:), se(3, :),'LineWidth', 2, 'cmap',green,'transparency',0.2,'alpha');
h = double(allPs<0.05); 
hb = h; hb(h==0) = nan; hb(hb==1) = -.015; 
set(gca, fontsize=25)
plot(times, hb,'k',  LineWidth=7)
set(gca, 'xlim', [-.25 1.7],'ylim', [-.02 .03]); 
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);


exportgraphics(gcf, ['myP.png'], 'Resolution',300)


%% 

c2u = 2; 
c1 = squeeze(allData_1(:, c2u, :)); 
c2 = squeeze(allData_2(:, c2u, :));

[h p ci ts] = ttest(c1, c2); 
%disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ';' ' p = ' num2str(p, 3)]);

plot(ts.tstat); hold on; 
plot(h)









%% EXTRACT STA LEVELS IN TMP or PFC CLUSTER and correlate with behavior


clear, clc

tyTR = 3; 
minTrN = 4;
remOutliers = 0; 
subtractionH = 0; % the difference between ACQ and EXT. If set to zero, only correlates with extinction (c1)

%tP = 19:31; %21:28 PFC effect %19:31 PFC effec ACQvsEXT
tP = 18:25; %TMP effect

paths = load_paths_EXT;

c1 = 'trlSTA_TMP_CTE_1-44_1_0_500-50';
c2 = 'trlSTA_TMP_CTA_1-44_1_0_500-50';
%c2 = 'none';
sub2exc = [37];
[cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
allRHO = zeros(50, 1); 


for subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 

    
    if ~isempty(rsa2T1)        
        % % % % % % take only CS+- items
        if tyTR ~= 0
            rsa2T1IDs = rsa2T1IDs(rsa2T1IDs(:, 6) ==tyTR,:); 
            rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) ==tyTR, :); 

            rsa2T2IDs = rsa2T2IDs(rsa2T2IDs(:, 6) ==tyTR,:); 
            rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) ==tyTR, :); 

        end

        clear trSTA
        for triali = 1:size(rsa2T1, 1)
            cTR1 = squeeze(mean(rsa2T1(triali, tP), 2));
            cTR2 = squeeze(mean(rsa2T2(triali, tP), 2));
            
            if subtractionH
                trSTA(triali, :) = cTR1 - cTR2;
            else
                trSTA(triali, :) = cTR1;     
            end
        end
        r2check{subji, :} = rsa2T1IDs(:, 7);
        ratings2u = rsa2T1IDs(:, 7);
        ratings2u2 = rsa2T2IDs(:, 7);
  
        % remove nans in both
        nanIds = isnan(ratings2u); 
        ratings2u(nanIds) = []; 
        trSTA(nanIds) = []; 
        

        if remOutliers
            %[B, tF1] = rmoutliers(trSTA, 'percentiles', [5 95]); 
            %[B, tF2] = rmoutliers(ratings2u, 'percentiles', [5 95]); 
            [B, tF1] = rmoutliers(trSTA, 'mean', ThresholdFactor=2); 
            [B, tF2] = rmoutliers(ratings2u, 'mean', ThresholdFactor=2); 
            tF3 = tF1 | tF2; 
            nOut(subji, :) = sum(tF3); 
            trSTA(tF3) = []; 
            ratings2u(tF3) = []; 
        end

        nTr(subji, :) = length(trSTA); 
        if ~isempty(trSTA)
            allRHO(subji, : ) = corr (trSTA, ratings2u, 'type', 's');
        end
       
    end
end

sub2exc2 = find(nTr < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc = unique(sub2exc3); 


allRHO(sub2exc) = []; 
nTr(allRHO==0 | isnan(allRHO)) = []; 
allRHO(allRHO==0 | isnan(allRHO)) = []; 
nTr(isnan(allRHO)) = []; 
allRHO(isnan(allRHO)) = []; 
allRHO = atanh(allRHO); 


[h p ci ts] = ttest(allRHO); 
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ';' ' p = ' num2str(p, 3)]);

figure(); set (gcf, 'Position', [50 50 300 500]);
initialColorOrder = get(gca,'ColorOrder'); 
scatter(1, allRHO, 75, initialColorOrder(1, :), 'filled'); hold on
scatter(1, allRHO, 200, initialColorOrder(2, :), 'ko', 'LineWidth', 1); hold on
plot(get(gca, 'xlim'), [0 0], 'k:', linewidth=2)
set(gca, 'ylim', [-1 1], 'xlim', [.5 1.5], 'xtick', [], 'xticklabel', [])
%title (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)])
set(gca, fontsize=24)


exportgraphics(gcf, 'myP.png', 'Resolution', 300);


%% EXTRACT STA LEVELS IN TMP or PFC CLUSTER and correlate with behavior for all three trial types


clear, clc


for trtyi = 1:3

    clearvars -except allRHO trtyi allRHOTR
    tyTR = trtyi; 
    minTrN = 4;

    %tP = 19:31; %21:28 PFC effect %19:31 PFC effec ACQvsEXT
    tP = 18:25; %TMP effect
    
    paths = load_paths_EXT;
    
    c1 = 'trlSTA_TMP_CTE_1-44_1_0_500-50';
    c2 = 'none';
    sub2exc = [37]';
    [cond1 cond2] = determine_conds_EXT(c1, c2, paths); 
    allRHO = zeros(50, 1); 
    
    
    for subji = 1:50
    
        rsa2T1 = cond1{subji, 1}; 
        rsa2T1IDs = cond1{subji, 2}; 
    
        
        if ~isempty(rsa2T1)        
            % % % % % % take only CS+- items
            if tyTR ~= 0
                rsa2T1IDs = rsa2T1IDs(rsa2T1IDs(:, 6) ==tyTR,:); 
                rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) ==tyTR, :); 
            end
    
            clear trSTA
            for triali = 1:size(rsa2T1, 1)
                cTR = squeeze(mean(rsa2T1(triali, tP), 2));
                trSTA(triali, :) = cTR;     
            end
            r2check{subji, :} = rsa2T1IDs(:, 7);
            ratings2u = rsa2T1IDs(:, 7);
    
            % remove nans in both
            nanIds = isnan(ratings2u); 
            ratings2u(nanIds) = []; 
            trSTA(nanIds) = []; 
            
    
            %[B, tF1] = rmoutliers(trSTA, 'percentiles', [5 95]); 
            %[B, tF2] = rmoutliers(ratings2u, 'percentiles', [5 95]); 
            [B, tF1] = rmoutliers(trSTA, 'mean', ThresholdFactor=2); 
            [B, tF2] = rmoutliers(ratings2u, 'mean', ThresholdFactor=2); 
            tF3 = tF1 | tF2; 
            nOut(subji, :) = sum(tF3); 
            trSTA(tF3) = []; 
            ratings2u(tF3) = []; 
    
            nTr(subji, :) = length(trSTA); 
            if ~isempty(trSTA)
                allRHO(subji, : ) = corr (trSTA, ratings2u, 'type', 's');
            end
           
        end
    end
    allRHOTR(:, trtyi) = allRHO; 

end

%% remove all zeros and nans

allRHOTR = allRHOTR(:, 1:2)
allRHOTR(any(isnan(allRHOTR), 2), :) = [];
allRHOTR = allRHOTR(any(allRHOTR,2),:);
allRHOTR(4,:) = [];

[h p ci ts] = ttest(allRHOTR(:, 1), allRHOTR(:, 2))


%%
colors = brewermap([100], '*Set1'); 

% rather than a square plot, make it thinner
violinPlot(allRHOTR, 'histOri', 'center', 'widthDiv', [2 1], 'showMM', 0, ...
    'color',  mat2cell(colors(1, : ), 1)); 
 
%set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'low', 'high'}, 'xlim', [0.2 1.8]);
ylabel('Value'); xlabel('Data');
xticks = get(gca, 'xtick');
for b = 1:3
    [~, pval] = ttest(allRHOTR(:, b));
    yval = max(allRHOTR(:, b)) * 1.2; % plot this on top of the bar
    
    mysigstar(gca, xticks(b), yval, pval);
    % if mysigstar gets just 1 xpos input, it will only plot stars
end


%% correlate context specificity with item stability across experimental phases Ext ->Test or ACQ > TEst


clear, clc


tyTR = 0; %0 = all trials
minTrN = 8; 


c1 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_TMP_CET_1-44_1_0_500-50'; 

printClust = 1; 
print1Clust = 0;


paths = load_paths_EXT;

sub2exc = [37]';
[cond1 cond2] = determine_conds_EXT(c1, c2, paths); 


nTimepoints = 51; 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

allRHO = zeros(50, bins, bins); 
nTrials = zeros(50, 1); 

f2s = [c1(1:10) '_' c2(1:10)]; 

parfor subji = 1:50

    rsa2T1 = cond1{subji, 1}; 
    rsa2T1IDs = cond1{subji, 2}; 

    rsa2T2 = cond2{subji, 1}; 
    rsa2T2IDs = cond2{subji, 2}; 
    
    if ~isempty(rsa2T1) & ~isempty(rsa2T2)
        [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
        rsa2T1 = rsa2T1(i1, :); 
        rsa2T1IDs = rsa2T1IDs(i1,:); 
        rsa2T2 = rsa2T2(i2, :); 
        rsa2T2IDs = rsa2T2IDs(i2,:); 


        % % % % % % take only CS+- items
        if tyTR ~= 0
            rsa2T1IDs = rsa2T1IDs(rsa2T1IDs(:, 6) ==tyTR,:); 
            rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) ==tyTR, :); 
            rsa2T2IDs = rsa2T2IDs(rsa2T2IDs(:, 6) ==tyTR,:); 
            rsa2T2 = rsa2T2(rsa2T2IDs(:, 6) ==tyTR, :); 
            % rsa2T1IDs = rsa2T1IDs(rsa2T1IDs(:, 6) ==2 | rsa2T1IDs(:, 6) ==3 ,:); 
            % rsa2T1 = rsa2T1(rsa2T1IDs(:, 6) ==2 | rsa2T1IDs(:, 6) ==3, :); 
            % rsa2T2IDs = rsa2T2IDs(rsa2T1IDs(:, 6) ==2 | rsa2T1IDs(:, 6) ==3,:); 
            % rsa2T2 = rsa2T2(rsa2T1IDs(:, 6) ==2 | rsa2T1IDs(:, 6) ==3, :); 
        end



        nTrials(subji, :) = size(rsa2T1, 1); 
    
        if ~isempty(rsa2T1)
             for timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 
        
                for timej = 1:bins
                    timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                    rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
           
                    allRHO(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
                end
             end
        end
    end
end

sub2exc2 = find(nTrials < minTrN); 
sub2exc3 = [sub2exc, sub2exc2']; 
sub2exc = unique(sub2exc3); 
allRHO(sub2exc, :, :) = []; 

[h p ci ts] = ttest(atanh(allRHO));
h = squeeze(h); t = squeeze(ts.tstat); 

%h(1:5, :) = 0;
%h(:, 1:5) = 0;

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs = allSTs(id); 
end

if ~printClust
    h = zeros(size(h)); 
end
if print1Clust
    h(clustinfo.PixelIdxList{id}) = 1; 
end

d2p = squeeze(mean(allRHO)); 
figure(); colormap (brewermap([100], '*Spectral'))
contourf(myresizem(t, 10), 50, 'linecolor', 'none'); axis square; hold on; 
contour( myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
plot([55 55],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(get(gca,'xlim'), [55 55],'k:', 'linewidth', 4); hold on; 
set(gca, 'xlim', [5 400], 'ylim', [5 400], 'FontSize', 24)
set(gca, 'xTick', [55 155 255 355], 'XTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'yTick', [55 155 255 355], 'YTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'clim', [-4 4])
ylabel('CTX_S_P_E')
exportgraphics(gcf, 'myP.png', 'Resolution', 300);




%% PERMUTATIONS
clearvars -except max_clust_obs cond1 cond2 sub2exc f2s
clc

% sub2exc inherited from previous code block

nPerm = 100;
t4P = 6:40; 

nTimepoints = length(t4P); 
win_width = 10; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );


paths = load_paths_EXT; 


tic 

for permi = 1:nPerm

    clearvars -except nTimepoints win_width mf bins cond1 cond2 permi nPerm max_clust_perm max_clust_obs t4P sub2exc f2s paths

    progress_in_console(permi)

    allRHOP = zeros(50, bins, bins); 
    parfor subji = 1:50
    
        rsa2T1 = cond1{subji, 1}; 
        rsa2T1IDs = cond1{subji, 2}; 
    
        rsa2T2 = cond2{subji, 1}; 
        rsa2T2IDs = cond2{subji, 2}; 
        
        if ~isempty(rsa2T1) & ~isempty(rsa2T2)
            [C i1 i2] = intersect(rsa2T1IDs(:, 1), rsa2T2IDs(:,1)); 
            rsa2T1 = rsa2T1(i1, t4P); 
            rsa2T1IDs = rsa2T1IDs(i1,:); 
    
            rsa2T2 = rsa2T2(i2, t4P); 
            rsa2T2IDs = rsa2T2IDs(i2,:); 

            ids4perm = randperm(size(rsa2T2, 1)); 
            rsa2T2 = rsa2T2(ids4perm, :); 
            rsa2T2IDs = rsa2T2IDs(ids4perm, :); 
        
            
             for timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                rsa2TT1 = mean(rsa2T1(:, timeBinsi), 2); 
    
                for timej = 1:bins
                    timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                    rsa2TT2 = mean(rsa2T2(:, timeBinsj), 2); 
           
                    allRHOP(subji, timei, timej) = corr(rsa2TT1, rsa2TT2, 'type', 's');
                end
            end
        end
    end
    
    allRHOP(sub2exc, :, :) = []; 
    
    [h p ci ts] = ttest(atanh(allRHOP));
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm(permi, : ) = allSTs(id); 
    else
        max_clust_perm(permi, : ) = 0; 
    end
    
    
end    

allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
%allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

f2s = [f2s '_' num2str(nPerm), 'p'];
save([paths.results.rsa_perm f2s], 'max_clust_perm', 'p', 'nPerm', 'max_clust_obs')


toc








%%


%allAb = max_clust_perm(abs(max_clust_perm) > abs(max_clust_obs));
allAb = max_clust_perm(max_clust_perm > max_clust_obs);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm




%% 

histogram (max_clust_perm)


%% 