%% Correlate amygdala power with different trial-level metrics
%% 
clear, clc



%timeAmyPow    = 6:15; 
%timeRSA     = 6-3:15-3; 

timeAmyPow  = 17:24;
timeRSA     = 3:12; 
%timeRSA     = 15:22; 



paths = load_paths_EXT;
%load ([paths.results.trial_based 'AMY_POW_3-8Hz'])
load ([paths.results.POWfromRT 'POW_AMY_C_100']); 
%load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_aHPC_CE_3-54_1_0_500-100'])
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

    amyPOW = POW{subji, 1}; 
    amyPOWIDs = double(string(POW{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 

    

    if ~isempty(amyPOW) & ~isempty(rsa2T)
        
        [i1 i2 i3] = intersect(amyPOWIDs(:, 1:9), rsa2TIDs(:,1:9), 'rows'); 
        amyPOW = squeeze(mean(mean(mean(amyPOW(i2, :, 3:8, timeAmyPow), 2), 3), 4)); 
        amyPOWIDs = amyPOWIDs(i2,:); 
        rsa2T = mean(rsa2T(i3, timeRSA), 2); 
        rsa2TIDs = rsa2TIDs(i3, :); 


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
hb = scatter([1], data.data, 50, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-1.7 1.7])
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'myP.png', 'Resolution', 300);




%% Correlate amygdala power IN CLUSTER with different trial-level metrics
%% 
clear, clc


%tp2use = 6:15; %18:23; 
tp2use = 18:23;

paths = load_paths_EXT;
load ([paths.results.trial_based 'AMY_POW_3-8Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
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
hb = scatter([1], data.data, 50, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-.7 .7])
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
load ([paths.results.trial_based 'AMY_POW_3-8Hz_TR'])
%load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])


%load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_aHPC_CE_3-54_1_0_500-100'])
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

    amyPOW1 = allPOWAMY{subji, 1}; 
    amyPOWIDs1 = double(string(allPOWAMY{subji, 2})); 

    rsa2T = cond2u{subji, 1}; 
    rsa2TIDs = cond2u{subji, 2}; 

    if ~isempty(amyPOW1) & ~isempty(rsa2T)
        
        for timei = 1:26
            [i1 i2] = intersect(amyPOWIDs1(:, 1), rsa2TIDs(:,1)); 
            amyPOW = amyPOW1(i2, :); 
            amyPOWIDs = amyPOWIDs1(i2,:); 

            rsa2TT = mean(rsa2T(:, timei), 2); 
    
            % % % z-score amygdala power separately for CS+ and CS-
            amyPOWCSp = amyPOW(amyPOWIDs(:, 8) == 1); 
            amyPOWCSm = amyPOW(amyPOWIDs(:, 8) == 0); 
            amyPOWCSp = (amyPOWCSp - mean(amyPOWCSp, 'omitnan')) ./ std(amyPOWCSp, 'omitnan');
            amyPOWCSm = (amyPOWCSm - mean(amyPOWCSm, 'omitnan')) ./ std(amyPOWCSm, 'omitnan');
            amyPOW(amyPOWIDs(:, 8) == 1) = amyPOWCSp; 
            amyPOW(amyPOWIDs(:, 8) == 0) = amyPOWCSm; 
            % % % z-score RSA2T metric separately for CS+ and CS-
            rsa2TCSp = rsa2TT(rsa2TIDs(:, 8) == 1); 
            rsa2TCSm = rsa2TT(rsa2TIDs(:, 8) == 0); 
            rsa2TCSp = (rsa2TCSp - mean(rsa2TCSp, 'omitnan')) ./ std(rsa2TCSp, 'omitnan');
            rsa2TCSm = (rsa2TCSm - mean(rsa2TCSm, 'omitnan')) ./ std(rsa2TCSm, 'omitnan');
            rsa2TT(rsa2TIDs(:, 8) == 1) = rsa2TCSp; 
            rsa2TT(rsa2TIDs(:, 8) == 0) = rsa2TCSm; 
    
            nanIds = isnan(amyPOW); 
            amyPOW(nanIds) = []; 
            rsa2TT(nanIds) = []; 
    
    
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
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

times = -.5:.1:2
figure()
plot(times, t, Linewidth=2); hold on; 
plot(times, hb, Linewidth=2)
set(gca, 'Fontsize', 16)
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);
exportgraphics(gcf, ['myP.png'], 'Resolution',150)



%% Correlate DIFFERENT TRIAL LEVEL METRICS
%% 
clear, clc

time1     = 14:20; 
time2     = 11:17; 

paths = load_paths_EXT;

load ([paths.results.trial_based 'trlSTA_AMY_CE_3-54_1_0_500-100'])
cond1 = itstaTRALL; 

%load ([paths.results.trial_based 'trlSTA_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_HPC_CE_3-54_1_0_500-100'])

load ([paths.results.trial_based 'trlSTA_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlSTA_OCC_CE_3-54_1_0_500-100'])
cond2 = itstaTRALL; 



%load ([paths.results.trial_based 'trlCTX_PFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_aHPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_HPC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_AMY_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_TMP_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OFC_CE_3-54_1_0_500-100'])
%load ([paths.results.trial_based 'trlCTX_OCC_CE_3-54_1_0_500-100'])
%cond2 = ctxTRALL; 


if length(cond2) < 50 
    cond2{50,2}= []; 
end



for subji = 1:50

    cond1S = cond1{subji, 1}; 
    cond1IDs = double(string(cond1{subji, 2})); 
    cond2S = cond2{subji, 1}; 
    cond2IDs = double(string(cond2{subji, 2})); 

    if ~isempty(cond1S) & ~isempty(cond2S)
        
        [i1 i2 i3]  = intersect(cond1IDs(:, 1:9), cond2IDs(:,1:9), 'rows'); 
        cond1P      = mean(cond1S(i2, time1), 2); 
        cond1PIDs   = cond1IDs(i2,:); 
        cond2P      = mean(cond2S(i3, time2), 2); 
        cond2PIDs   = cond2IDs(i3, :); 

        
        % % % z-score RSA2T metric separately for CS+ and CS-
        cond1CSp = cond1P(cond1PIDs(:, 8) == 1); 
        cond1CSm = cond1P(cond1PIDs(:, 8) == 0); 
        cond1CSp = (cond1CSp - mean(cond1CSp, 'omitnan')) ./ std(cond1CSp, 'omitnan');
        cond1CSm = (cond1CSm - mean(cond1CSm, 'omitnan')) ./ std(cond1CSm, 'omitnan');
        cond1P(cond1PIDs(:, 8) == 1) = cond1CSp; 
        cond1P(cond1PIDs(:, 8) == 0) = cond1CSm; 

        cond2CSp = cond2P(cond2PIDs(:, 8) == 1); 
        cond2CSm = cond2P(cond2PIDs(:, 8) == 0); 
        cond2CSp = (cond2CSp - mean(cond2CSp, 'omitnan')) ./ std(cond2CSp, 'omitnan');
        cond2CSm = (cond2CSm - mean(cond2CSm, 'omitnan')) ./ std(cond2CSm, 'omitnan');
        cond2P(cond2PIDs(:, 8) == 1) = cond2CSp; 
        cond2P(cond2PIDs(:, 8) == 0) = cond2CSm; 

        % nanIds = isnan(amyPOW); 
        % amyPOW(nanIds) = []; 
        % rsa2T(nanIds) = []; 


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


        if length(cond1P) > 5 % at least five trials for the correlation analysis
            allRHO(subji, :) = corr(cond1P, cond2P, 'type', 's');
        end

    else
        
        allRHO(subji, :) = nan; 

    end

end


allRHO(isnan(allRHO)) = []; 
[h p ci ts] = ttest(atanh(allRHO));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


%%plot one bar
clear data
data.data = atanh(allRHO); 

figure(); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 50, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 25, 'linew',1, 'xlim', [0 2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2);
set(gca, 'LineWidth', 2);
set(gca, 'ylim', [-.8 .8])
box on; 
[h p ci ts] = ttest (data.data);
%res2title = ['t = ' num2str(t.tstat, '%.3f') '  ' ' p = ' num2str(p, '%.3f')]; 
res2title = ['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)];
disp (res2title);

%title(res2title)

exportgraphics(gcf, 'myP.png', 'Resolution', 300);












































%% 