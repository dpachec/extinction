%% Correlate amygdala power IN CLUSTER with different trial-level metrics
%% 
clear, clc

%tp2use = 6:15; %18:23; 
tp2use = 18:23;

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_3-54Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])


load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_3-54_1_0_500-100'])
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
        rsa2TCSp = rsa2T(rsa2TIDs(:, 8) == 1); 
        rsa2TCSm = rsa2T(rsa2TIDs(:, 8) == 0); 
        rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
        rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
        rsa2T(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
        rsa2T(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 

        nanIds = isnan(amyPOW); 
        amyPOW(nanIds) = []; 
        rsa2T(nanIds) = []; 


        %[B, tF1] = rmoutliers(rsa2T, 'percentiles', [10 90]); 
        %[B, tF2] = rmoutliers(amyPOW, 'percentiles', [10 90]); 
        % [B, tF1] = rmoutliers(rsa2T, 'mean'); 
        % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
        % tF3 = tF1 | tF2; 
        % nOut(subji, :) = sum(tF3); 
        % rsa2T(tF3) = []; 
        % amyPOW(tF3) = []; 


        % % figure()
        % % scatter(amyPOW, rsa2T, 150, 'filled');
        % % h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
        % % C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
        % % allSlopes(subji, :) = C(2);
        % % allIntercepts(subji, :) = C(1);
        % % %set(gca, 'ylim', [-.15 .15], 'xlim', [-3 3], 'Fontsize', 24)
        % % set(gca, 'ylim', [-3 3], 'xlim', [-3 3], 'Fontsize', 24)


        if length(amyPOW) > 5 % at least five trials for the correlation analysis
            allRHO(subji, :) = corr(amyPOW, rsa2T, 'type', 's');
        end

    else
        
        allRHO(subji, :) = nan; 

    end

end


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
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-1 1])
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'myP.png', 'Resolution', 300);





%% Correlate amygdala power IN CLUSTER with different trial-level metrics (ALL TIME POINTS)
clear, clc



paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_3-54Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_9-54_1_0_500-100'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_PFC_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_9-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_OCC_CE_9-54_1_0_500-100'])




if exist('itstaTRALL')
    cond2u = itstaTRALL; 
end
if exist('ctxTRALL')
    cond2u = ctxTRALL; 
end

if length(cond2u) < 50 
    cond2u{50,2}= []; 
end


nTimepoints = 26; %%all2all file is stored within this cell array
win_width = 5; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

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
    
            
    
    
            %[B, tF1] = rmoutliers(rsa2T, 'percentiles', [10 90]); 
            %[B, tF2] = rmoutliers(amyPOW, 'percentiles', [10 90]); 
            % [B, tF1] = rmoutliers(rsa2T, 'mean'); 
            % [B, tF2] = rmoutliers(amyPOW, 'mean'); 
            % tF3 = tF1 | tF2; 
            % nOut(subji, :) = sum(tF3); 
            % rsa2T(tF3) = []; 
            % amyPOW(tF3) = []; 
    
    
            % % figure()
            % % scatter(amyPOW, rsa2T, 150, 'filled');
            % % h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
            % % C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
            % % allSlopes(subji, :) = C(2);
            % % allIntercepts(subji, :) = C(1);
            % % %set(gca, 'ylim', [-.15 .15], 'xlim', [-3 3], 'Fontsize', 24)
            % % set(gca, 'ylim', [-3 3], 'xlim', [-3 3], 'Fontsize', 24)
    
    
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
    max_clust_obs = allSTs(id); 
end


%h(1:10) = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

mAR = mean(allRHO);
stdAR = std(allRHO); 
seAR = stdAR / sqrt(size(allRHO, 1))

times = -.25:.1:1.9
%times = -.25:.1:2.2
figure()

%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'Set3')*0.75;
shadedErrorBar(times, mAR, seAR, {'Color',colors2use(1,:)}, 1); hold on; 

%set(gca, 'ylim', [-.25 .1], 'xlim', [-.25 1.8]) % for AMY HPC STABILITY
set(gca, 'ylim', [-.2 .2], 'xlim', [-.25 1.8]) % 

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








%% Correlate different trial-level metrics at ALL TIME POINTS
clear, clc

printClust = 0; 
print1Clust = 1; 

paths = load_paths_EXT;

cond1 = load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%cond2 = load([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])


%cond2 = load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlSTA_HPC_CE_3-54_1_0_500-100'])
cond2 = load ([paths.results.trial_based 'trlSTA_PFC_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlSTA_TMP_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlSTA_OFC_CE_3-54_1_0_500-100'])
%cond2 = load ([paths.results.trial_based 'trlSTA_OCC_CE_3-54_1_0_500-100'])



if isfield(cond1, 'ctxTRALL')
    cond1 = cond1.ctxTRALL; 
else
    cond1 = cond1.itstaTRALL; 
end

if isfield(cond2, 'ctxTRALL')
    cond2 = cond2.ctxTRALL; 
else
    cond2 = cond2.itstaTRALL; 
end

if length(cond1) < 50 
    cond1{50,2}= []; 
end
if length(cond2) < 50 
    cond2{50,2}= []; 
end


nTimepoints = 26; 
win_width = 5; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

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

% % % % only for tmp
%sub2exc = [21];
%allRHO(sub2exc, :, :) = []; 

ids2K = any(allRHO,2); 
ids2K = ids2K(:, 1); 
allRHO = allRHO(ids2K,:,:);
[h p ci ts] = ttest(atanh(allRHO));
h = squeeze(h); t = squeeze(ts.tstat); 



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
contourf(myresizem(t, 10), 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
plot([25 25],get(gca,'ylim'),'k:', 'linewidth', 4); hold on; 
plot(get(gca,'xlim'), [25 25],'k:', 'linewidth', 4); hold on; 
set(gca, 'xlim', [1 200], 'ylim', [1 200], 'FontSize', 30)
set(gca, 'xTick', [25 75 125 175], 'XTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'yTick', [25 75 125 175], 'YTickLabel', {'0', '0.5', '1', '1.5'})
set(gca, 'clim', [-4 4])

exportgraphics(gcf, 'myP.png', 'Resolution', 300);



















%% PERMUTATIONS
clearvars -except max_clust_obs cond1 cond2
clc

nPerm = 1000;
t4P = 3:20; 

nTimepoints = 18; 
win_width = 5; 
mf = 1; 
bins =  floor ( (nTimepoints/mf)- win_width/mf+1 );

tic

for permi = 1:nPerm

    clearvars -except nTimepoints win_width mf bins cond1 cond2 permi nPerm max_clust_perm max_clust_obs t4P
    for subji = 1:50
    
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
        
            
             parfor timei = 1:bins 
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
    
    % % % % only for tmp
    %sub2exc = [21];
    %allRHO(sub2exc, :, :) = []; 
    
    ids2K = any(allRHO,2); 
    ids2K = ids2K(:, 1); 
    allRHO = allRHO(ids2K,:,:); % limit also time 
    [h p ci ts] = ttest(atanh(allRHO));
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
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm




toc


%%










































%% 