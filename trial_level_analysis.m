%% Correlate amygdala power IN CLUSTER with different trial-level metrics
%% 
clear, clc

%tp2use = 36:46;
%tp2use = 20:32; % check OFC
tp2use = 33:41; % > AMY theta power effect


paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_1-44Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_1-44_1_0_500-50'])
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
        

        [i1 i2] = intersect(amyPOWIDs(:, 1), rsa2TIDs(:,1)); 
        amyPOW = amyPOW(i2, :); 
        amyPOWIDs = amyPOWIDs(i2,:); 
        rsa2T = mean(rsa2T(:, tp2use), 2); 


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


        if length(amyPOW) > 5 % at least five trials for the correlation analysis
            allRHO(subji, :) = corr(amyPOW, rsa2T, 'type', 's');
        end

    else
        
        allRHO(subji, :) = nan; 

    end

end

%allRHO(sub2exc) = nan; 
%allRHO(isnan(allRHO)) = []; 
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





%% Correlate amygdala POWER IN CLUSTER with different trial-level metrics (ALL TIME POINTS)
clear, clc

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

sub2exc = [37]; 
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
















%% Correlate different trial-level metrics at ALL TIME POINTS NEW TIMES
clear, clc


minTrN = 8; 

printClust = 0; 
print1Clust = 0; 

paths = load_paths_EXT;

c1 = 'trlSTA_PFC_CE_1-44_1_0_500-50'; 
c2 = 'trlSTA_AMY_CE_1-44_1_0_500-50'; 

%c1 = 'trlCTX_AMY_CE_1-44_1_0_500-50'; 
%c2 = 'trlCTX_PFC_CE_1-44_1_0_500-50'; 

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

        
        % % % z-score both trial-level metrics separately for CS+ and CS-
        rsa2TCSp = rsa2T1(rsa2T1IDs(:, 8) == 1,:); 
        rsa2TCSm = rsa2T1(rsa2T1IDs(:, 8) == 0, :); 
        rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        rsa2T1(rsa2T1IDs(:, 8) == 1,:) = rsa2TCSp; 
        rsa2T1(rsa2T1IDs(:, 8) == 0,:) = rsa2TCSm; 


        rsa2TCSp = rsa2T2(rsa2T2IDs(:, 8) == 1,:); 
        rsa2TCSm = rsa2T2(rsa2T2IDs(:, 8) == 0,:); 
        rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        rsa2T2(rsa2T2IDs(:, 8) == 1,:) = rsa2TCSp; 
        rsa2T2(rsa2T2IDs(:, 8) == 0,:) = rsa2TCSm; 


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

% h = zeros(size(h)); 
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







%% correlate CTX SPE in PFC across SUBJECTS with REINST in AMY AND TMP (across subjects)


clear, clc

trltype = 2;
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

rsa2TT4 = rsa2TT2 - rsa2TT3; 
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







%% 

histogram(max_clust_perm); hold on; 
scatter(max_clust_obs,0, 5000, 'r.')




%% EXTRACT STA LEVELS IN TMP or PFC CLUSTER and correlate with behavior


clear, clc

tyTR = 2; 
minTrN = 8;
remOutliers = 0; 
subtractionH = 1; % the difference between ACQ and EXT. If set to zero, only correlates with extinction (c1)

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
%nTr(allRHO==0 | isnan(allRHO)) = []; 
%allRHO(allRHO==0 | isnan(allRHO)) = []; 
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