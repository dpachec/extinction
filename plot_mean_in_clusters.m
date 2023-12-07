%% 
%% plot mean in cluster 
% % % % load cluster 

clear 
load HPC_SCA-DCA_px1

p2u = 1; 
paths = load_paths_EXT; 
f2sav =  'POW_HPC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond1a cond1b] = rem_half_matrix(cond1, cond2);

f2sav =  'POW_HPC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond2a cond2b] = rem_half_matrix(cond1, cond2);


for subji = 1:size(cond1, 1)

    c1a = squeeze(cond1a(subji, :, :)); 
    c1b = squeeze(cond1b(subji, :, :)); 
    mc1(subji, :) = mean(c1a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c1b(clustinfo.PixelIdxList{p2u}), 'all'); 
    c2a = squeeze(cond2a(subji, :, :)); 
    c2b = squeeze(cond2b(subji, :, :)); 
    mc2(subji, :) = mean(c2a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c2b(clustinfo.PixelIdxList{p2u}), 'all'); 
    

end


data.data = [mc1 mc2];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.1 .05] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

%[h p ci t] = ttest (data.data(:,1), data.data(:,2));
[h p ci t] = ttest (data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)
%close all;   


%% plot mean in cluster 
% % % % load cluster 

clear 
load HPC_SCE-DCE_px1

p2u = 1; 
paths = load_paths_EXT; 
f2sav =  'POW_HPC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond1a cond1b] = rem_half_matrix(cond1, cond2);

f2sav =  'POW_HPC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond2a cond2b] = rem_half_matrix(cond1, cond2);


for subji = 1:size(cond1, 1)

    c1a = squeeze(cond1a(subji, :, :)); 
    c1b = squeeze(cond1b(subji, :, :)); 
    mc1(subji, :) = mean(c1a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c1b(clustinfo.PixelIdxList{p2u}), 'all'); 
    c2a = squeeze(cond2a(subji, :, :)); 
    c2b = squeeze(cond2b(subji, :, :)); 
    mc2(subji, :) = mean(c2a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c2b(clustinfo.PixelIdxList{p2u}), 'all'); 
    

end


data.data = [mc1 mc2];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.05 .1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

%[h p ci t] = ttest (data.data(:,1), data.data(:,2));
[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)
%close all;   


%% plot mean in cluster 
% % % % load cluster 

clear 
load PFC_SCA-DCA_px3

p2u = 1; 
paths = load_paths_EXT; 
f2sav =  'POW_PFC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond1a cond1b] = rem_half_matrix(cond1, cond2);

f2sav =  'POW_PFC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond2a cond2b] = rem_half_matrix(cond1, cond2);


for subji = 1:size(cond1, 1)

    c1a = squeeze(cond1a(subji, :, :)); 
    c1b = squeeze(cond1b(subji, :, :)); 
    mc1(subji, :) = mean(c1a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c1b(clustinfo.PixelIdxList{p2u}), 'all'); 
    c2a = squeeze(cond2a(subji, :, :)); 
    c2b = squeeze(cond2b(subji, :, :)); 
    mc2(subji, :) = mean(c2a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c2b(clustinfo.PixelIdxList{p2u}), 'all'); 
    

end


data.data = [mc1 mc2];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.1 .2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
%[h p ci t] = ttest (data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)
%close all;   


%% plot mean in cluster 
% % % % load cluster 

clear 
load PFC_SCE-DCE_px2

p2u = 1; 
paths = load_paths_EXT; 
f2sav =  'POW_PFC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond1a cond1b] = rem_half_matrix(cond1, cond2);

f2sav =  'POW_PFC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
load ([ paths.results.rsa f2sav '.mat']);
ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
[cond2a cond2b] = rem_half_matrix(cond1, cond2);


for subji = 1:size(cond1, 1)

    c1a = squeeze(cond1a(subji, :, :)); 
    c1b = squeeze(cond1b(subji, :, :)); 
    mc1(subji, :) = mean(c1a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c1b(clustinfo.PixelIdxList{p2u}), 'all'); 
    c2a = squeeze(cond2a(subji, :, :)); 
    c2b = squeeze(cond2b(subji, :, :)); 
    mc2(subji, :) = mean(c2a(clustinfo.PixelIdxList{p2u}), 'all') - mean(c2b(clustinfo.PixelIdxList{p2u}), 'all'); 
    

end


data.data = [mc1 mc2];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.1 .2] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
%[h p ci t] = ttest (data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, 'myPNG.png', 'Resolution',150)
%close all;   


%% 