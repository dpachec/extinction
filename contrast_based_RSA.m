%% RSA IN LOOP 
%%
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast

clear , clc

listF2sav = {

'RSA_TMP_C_1-44_1_0_500-50_1_T_SICSPE-SICSME';
'RSA_PFC_C_1-44_1_0_500-50_1_T_SCE-DCE'; 

};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
        
    f2sav       = listF2sav{listi}; 
    cfg         = getParams_EXT(f2sav);
    
    
    ALLPOW = load ([paths.results.POWfromRT cfg.powF2load]); 
    
    
    for subji = 1:length(ALLPOW.POW)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        cfg.oneListPow = ALLPOW.POW{subji, 1};
        cfg.oneListIds = ALLPOW.POW{subji, 2};
        
        
        if ~isempty(cfg.oneListPow)

            out_contrasts = create_contrasts_EXT(cfg);
            %out_rsa(subji, :, :, :,:) = rsa_EXT(out_contrasts, cfg);
            out_rsa{subji, :} = rsa_EXT(out_contrasts, cfg);
            ids{subji,:} = out_contrasts.allIDs;

        end
    
        
    end

    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end





%% plot 2 lines from TG DIFFERENT TIMES TRIALS
clear

paths = load_paths_EXT; 
f2sav = 'RSA_OFC_C_1-44_1_0_500-50_1_T_SICSPE-SICSME';
%f2sav  = 'RSA_PFC_C_1-44_1_0_500-50_1_T_SCA-DCA';

load ([ paths.results.rsa f2sav '.mat']);

cfg = getParams_EXT(f2sav);
if length(cfg.contr2sav{1}) > 3 & strcmp(cfg.contr2sav{1}(1:5), 'SICSP') %ITEM STABILITY
    if strcmp(cfg.contr2sav{1}(6), 'E') 
        if strcmp(cfg.roi, 'TMP') 
            sub2exc = [19 25 31 27 37]; 
        elseif strcmp(cfg.roi, 'AMY') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'HPC') 
            sub2exc = [27 37];             
        elseif strcmp(cfg.roi, 'PFC') 
            sub2exc = [27 37];             
        elseif strcmp(cfg.roi, 'OFC') 
            sub2exc = [27 37];        
        elseif strcmp(cfg.roi, 'OCC') 
            sub2exc = [27 37];               
        end
    elseif strcmp(cfg.contr2sav{1}(6), 'A') 
        disp('Acquisition')
        sub2exc = [24 25 27 37]; 
    end
elseif strcmp(cfg.contr2sav{1}(1:2), 'SC') %CONTEXT
    if strcmp(cfg.contr2sav{1}(3), 'A') 
        if strcmp(cfg.roi, 'TMP') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'PFC') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'HPC') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'OCC') 
            sub2exc = [27 37];                  
        end
    elseif strcmp(cfg.contr2sav{1}(3), 'E')

    end
end

nTrials = compute_trial_number_EXT(ids); 

[out_rsa, ids2rem] = rem_nan_subj_EXT(out_rsa, sub2exc); 

cond1 = cellfun(@mean, cellfun(@(x) x{1}, out_rsa, 'un', 0), 'un', 0); 
cond2 = cellfun(@mean, cellfun(@(x) x{2}, out_rsa, 'un', 0), 'un', 0); 

cond1 = cat(1, cond1{:}); 
cond2 = cat(1, cond2{:}); 


diff = cond1-cond2; 

for subji = 1:size(cond1, 1)
   cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
   cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
end       
cond1 = cond1B(:, 1:40); cond2 = cond2B(:, 1:40); 



d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1));
se2 = d2pstd2/sqrt(size(cond1, 1));

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
else
    tObs = 0; 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
%h = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 


%times = (-.5:.01:2) + .25;
times = -.2:.05:1.75;
% to check the times 
h(2,:) = times; 
figure(); 
%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'*Accent')*0.75;
%colors2use = brewermap([6],'Accent')*0.75;

shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(2,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(1,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)
%set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 28);
set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 22);
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);

%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%
%title(f2sav, 'Interpreter','none')
%exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)








%% Permutations (2D) TRIAL BASED

nPerm = 100;
    
clearvars -except out_rsa nPerm tObs
for permi = 1:nPerm
    
    permi
    out_rsa_perm = randomize_trials_EXT(out_rsa); 
    cond1P = cellfun(@mean, cellfun(@(x) x{1}, out_rsa_perm, 'un', 0), 'un', 0); 
    cond2P = cellfun(@mean, cellfun(@(x) x{2}, out_rsa_perm, 'un', 0), 'un', 0);     
    cond1P = cat(1, cond1P{:}); 
    cond2P = cat(1, cond2P{:}); 
    for subji = 1:size(cond1P, 1)
       cond1BP(subji, :) = diag(squeeze(cond1P(subji, :, :)));
       cond2BP(subji, :) = diag(squeeze(cond2P(subji, :, :)));
    end       
    cond1P = mean(cond1BP); cond2P = mean(cond2BP); 
    [h p ci ts] = ttest(cond1P, cond2P); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
nPerm = length(max_clust_sum_perm); 
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% plot 2 lines from TG DIFFERENT TIMES (ONLY WITH THE AVERAGES)
clear

paths = load_paths_EXT; 
%f2sav = 'RSA_OFC_C_1-44_1_0_500-50_1_SICSPE-SICSME';
f2sav  = 'RSA_PFC_C_1-44_1_0_500-50_1_SCA-DCA';

load ([ paths.results.rsa f2sav '.mat']);

cfg = getParams_EXT(f2sav);
if length(cfg.contr2sav{1}) > 3 & strcmp(cfg.contr2sav{1}(1:5), 'SICSP') %ITEM STABILITY
    if strcmp(cfg.contr2sav{1}(6), 'E') 
        if strcmp(cfg.roi, 'TMP') 
            sub2exc = [19 25 31 27 37]; 
        elseif strcmp(cfg.roi, 'AMY') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'HPC') 
            sub2exc = [27 37];             
        elseif strcmp(cfg.roi, 'PFC') 
            sub2exc = [27 37];             
        elseif strcmp(cfg.roi, 'OFC') 
            sub2exc = [27 37];        
        elseif strcmp(cfg.roi, 'OCC') 
            sub2exc = [27 37];               
        end
    elseif strcmp(cfg.contr2sav{1}(6), 'A') 
        disp('Acquisition')
        sub2exc = [24 25 27 37]; 
    end
elseif strcmp(cfg.contr2sav{1}(1:2), 'SC') %CONTEXT
    if strcmp(cfg.contr2sav{1}(3), 'A') 
        if strcmp(cfg.roi, 'TMP') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'PFC') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'HPC') 
            sub2exc = [27 37]; 
        elseif strcmp(cfg.roi, 'OCC') 
            sub2exc = [27 37];                  
        end
    elseif strcmp(cfg.contr2sav{1}(3), 'E')

    end
end

nTrials = compute_trial_number_EXT(ids); 

[out_rsa, ids2rem] = rem_nan_subj_EXT(out_rsa, sub2exc); 

cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 



cond1(ids2rem, :, :) = []; 
cond2(ids2rem, :, :) = []; 

diff = cond1-cond2; 

for subji = 1:size(cond1, 1)
   cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
   cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
end       
cond1 = cond1B(:, 1:40); cond2 = cond2B(:, 1:40); 



d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1));
se2 = d2pstd2/sqrt(size(cond1, 1));

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
else
    tObs = 0; 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
%h = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 


%times = (-.5:.01:2) + .25;
times = -.2:.05:1.75;
% to check the times 
h(2,:) = times; 
figure(); 
%colors2use = brewermap([6],'*Set1')*0.75;
colors2use = brewermap([6],'*Accent')*0.75;
%colors2use = brewermap([6],'Accent')*0.75;

shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(2,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(1,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)
%set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 28);
set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 22);
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);

%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%
%title(f2sav, 'Interpreter','none')
%exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)










%% Permutations (2D)

nPerm = 1000;

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';


junts = cat(1, cond1(:, 6:40), cond2(:, 6:40));
%junts = cat(1, cond1(:, 3:20), cond2(:, 3:20));
%junts = cat(1, cond1(:, 3:22), cond2(:, 3:22));

[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);
    
clear max_clust_sum_perm h
for permi = 1:nPerm
    
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);

    cond1P = junts(fakeCondMapping == 0, :);
    cond2P = junts(fakeCondMapping == 1, :);

    diffC = cond1P - cond2P; 
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm



%% 
histogram (max_clust_sum_perm)



%% plot 2 lines from TG 100ms
clear


sub2exc = []; %

paths = load_paths_EXT; 
f2sav = 'RSA_TMP_C_3-54_1_0_500-100_1_SICSPE-SICSME';
%f2sav = 'RSA_TMP_C_1-44_1_0_500-100_1_SICSPE-SICSME';
%f2sav = 'RSA_PFC_C_1-44_1_0_500-100_1_SCE-DCE';

load ([ paths.results.rsa f2sav '.mat']);


cfg = getParams_EXT(f2sav);
if length(cfg.contr2sav{1}) > 3 & strcmp(cfg.contr2sav{1}(1:5), 'SICSP') 
    sub2exc = [19 25 31 27 37]; 
elseif strcmp(cfg.contr2sav{1}(1:2), 'SC')
    sub2exc = [27 37]; 
end


nTrials = compute_trial_number_EXT(ids); 

[out_rsa, ids2rem] = rem_nan_subj_EXT(out_rsa, sub2exc); 

cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids2rem, :, :) = []; 
cond2(ids2rem, :, :) = []; 

%cond1 = cond1 - mean(cond1(:, 1:2), 2); %subtract mean baseline
%cond2 = cond2 - mean(cond2(:, 1:2), 2); 

diff = cond1-cond2; 

for subji = 1:size(cond1, 1)
   cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
   cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
end       
cond1 = cond1B; cond2 = cond2B; 

% cond1(sub2exc, :) = []; 
% cond2(sub2exc, :) = []; 

d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1));
se2 = d2pstd2/sqrt(size(cond1, 1));

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
else
    tObs = 0; 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
%h = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

%times = (-.5:.01:2) + .25;
times = -.25:.1:2.3;
figure(); 
%colors2use = brewermap([6],'*Set1')*0.75;
%colors2use = brewermap([6],'*Accent')*0.75;
colors2use = brewermap([6],'Accent')*0.75;

shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(2,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(1,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)
%set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 28);
set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .04], 'Fontsize', 28);
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 3);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 3);

%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%
%title(f2sav, 'Interpreter','none')
%exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)




%% plot TG
clear, clc
paths = load_paths_EXT; 
 
f2sav = 'RSA_TMP_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
sub2exc = [];


load ([ paths.results.rsa f2sav '.mat']);

out_rsa(sub2exc, :, :,:) = []; 
ids = rem_nan_subj_EXT(out_rsa); 

%cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
%cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
cond1 = squeeze(out_rsa(:, 1, 1:20, 1:20)); 
cond2 = squeeze(out_rsa(:, 2, 1:20, 1:20)); 


cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
diff = cond1-cond2; 

[cond1 cond2] = rem_half_matrix(cond1, cond2);

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
%h(1:251,1:26) = 0;  % % % no clusters before baseline

clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    %[max2u id] = max((allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


clim = [-.03 .03];
plot_TG_map(m1, m2, h, t, f2sav, clim)
%exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)


%% Permutations (3D)

nPerm = 1000;

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';


junts = cat(1, cond1(:, 6:40, 6:40), cond2(:, 6:40, 6:40));
%junts = cat(1, cond1(:, 3:20, 3:20), cond2(:, 3:20, 3:20));
%junts = cat(1, cond1(:, 3:22, 3:22), cond2(:, 3:22, 3:22));

[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);
    
clear max_clust_sum_perm
for permi = 1:nPerm
    
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);

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
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm





%% CHECK CSP vs CSM during ACQ     VS    CSP vs CSM during EXT 2D
clear, clc
paths = load_paths_EXT; 

myR = 'TMP';

f2sav =   ['RSA_' myR '_C_1-44_1_0_500-50_1_SICSPA-SICSMA'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 
nTrialsA = compute_trial_number_EXT(ids); 
f2sav =   ['RSA_' myR '_C_1-44_1_0_500-50_1_SICSPE-SICSME'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 
nTrialsE = compute_trial_number_EXT(ids); 


cfg = getParams_EXT(f2sav);
if length(cfg.contr2sav{1}) > 3 & strcmp(cfg.contr2sav{1}(1:5), 'SICSP') 
    sub2exc = [19 24 25 31 27 37]; 
elseif strcmp(cfg.contr2sav{1}(1:2), 'SC')
    sub2exc = [27 37]; 
end





[out_rsa_ACQ, ids2remA] = rem_nan_subj_EXT(out_rsa_ACQ, sub2exc); 
cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); 
cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); 
cond1A(ids2remA, :, :) = []; 
cond2A(ids2remA, :, :) = []; 
diffA = cond1A-cond2A; 
diffA(sub2exc, :,:) = []; 

[out_rsa_EXT, ids2remE] = rem_nan_subj_EXT(out_rsa_EXT, sub2exc); 
cond1E = squeeze(out_rsa_EXT(:, 1, :, :)); 
cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); 
cond1E(ids2remE, :, :) = []; 
cond2E(ids2remE, :, :) = []; 
diffE = cond1E-cond2E;
diffE(sub2exc, :,:) = []; 

for subji = 1:size(diffA, 1)
   diffAB(subji, :) = diag(squeeze(diffA(subji, :, :)));
   diffEB(subji, :) = diag(squeeze(diffE(subji, :, :)));
end       
diffA = diffAB; diffE = diffEB; 



d2pm1	= squeeze(mean(diffA,'omitnan'));
d2pm2	= squeeze(mean(diffE,'omitnan'));
d2pstd1	= std(diffA, 'omitnan');
d2pstd2	= std(diffE, 'omitnan');
se1 = d2pstd1/sqrt(size(diffA, 1));
se2 = d2pstd2/sqrt(size(diffE, 1));


[h p ci ts] = ttest(diffA, diffE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
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
h(1:6) =0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.025; 

%times = (-.5:.01:2) + .25;
times = -.25:.05:2.25;
figure(); 
colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)

set(gca, 'xlim', [-.25 1.75],'ylim', [-.03 .045], 'Fontsize', 22);%
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);
%title(f2sav, 'Interpreter','none')
exportgraphics(gcf, ['myP.png'], 'Resolution',300)



%% CHECK CONTEXT DURING ACQ AND EXT > GENERALIZATION
clear, clc
paths = load_paths_EXT; 

myR = 'PFC';

f2sav =   ['RSA_' myR '_C_1-44_1_0_500-50_1_SCA-DCA'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 
nTrialsA = compute_trial_number_EXT(ids); 
f2sav =   ['RSA_' myR '_C_1-44_1_0_500-50_1_SCE-DCE'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 
nTrialsE = compute_trial_number_EXT(ids); 


cfg = getParams_EXT(f2sav);
if length(cfg.contr2sav{1}) > 3 & strcmp(cfg.contr2sav{1}(1:5), 'SICSP') 
    sub2exc = [19 24 25 31 27 37]; 
elseif strcmp(cfg.contr2sav{1}(1:2), 'SC')
    sub2exc = [27 37]; 
end

[out_rsa_ACQ, ids2remA] = rem_nan_subj_EXT(out_rsa_ACQ, sub2exc); 
cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); 
cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); 
cond1A(ids2remA, :, :) = []; 
cond2A(ids2remA, :, :) = []; 
diffA = cond1A-cond2A; 

[out_rsa_EXT, ids2remE] = rem_nan_subj_EXT(out_rsa_EXT, sub2exc); 
cond1E = squeeze(out_rsa_EXT(:, 1, :, :)); 
cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); 
cond1E(ids2remE, :, :) = []; 
cond2E(ids2remE, :, :) = []; 
diffE = cond1E-cond2E;

for subji = 1:size(diffA, 1)
   diffAB(subji, :) = diag(squeeze(diffA(subji, :, :)));
   diffEB(subji, :) = diag(squeeze(diffE(subji, :, :)));
end       
diffA = diffAB; diffE = diffEB; 



d2pm1	= squeeze(mean(diffA,'omitnan'));
d2pm2	= squeeze(mean(diffE,'omitnan'));
d2pstd1	= std(diffA, 'omitnan');
d2pstd2	= std(diffE, 'omitnan');
se1 = d2pstd1/sqrt(size(diffA, 1));
se2 = d2pstd2/sqrt(size(diffE, 1));


[h p ci ts] = ttest(diffA, diffE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
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
h(1:6) =0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.012; 

%times = (-.5:.01:2) + .25;
times = -.25:.05:2.25;
figure(); 
colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)

set(gca, 'xlim', [-.25 1.75],'ylim', [-.015 .02], 'Fontsize', 22);%
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);
%title(f2sav, 'Interpreter','none')
exportgraphics(gcf, ['myP.png'], 'Resolution',300)



%% CHECK CONTEXT DURING ACQ AND EXT > GENERALIZATION 2D
clear, clc
paths = load_paths_EXT; 

myR = 'PFC';

%f2sav =   ['RSA_' myR '_C_3-54_1_0_500-100_1_SICSPA-SICSMA'];
f2sav =   ['RSA_' myR '_C_3-54_1_0_500-100_1_SCA-DCA'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVA-DIDVA'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 

%f2sav =   ['RSA_' myR '_C_3-54_1_0_500-100_1_SICSPE-SICSME'];
f2sav =   ['RSA_' myR '_C_3-54_1_0_500-100_1_SCE-DCE'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVE-DIDVE'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 



% % % % %  remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)
    %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
    %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end
end

% cond1A = squeeze(out_rsa_ACQ(:, 1, 1:13, 1:13)); cond1A(ids, :, :) = []; 
% cond2A = squeeze(out_rsa_ACQ(:, 2, 1:13, 1:13)); cond2A(ids, :,:) = []; 
% cond1E = squeeze(out_rsa_EXT(:, 1, 1:13, 1:13)); cond1E(ids, :, :) = []; 
% cond2E = squeeze(out_rsa_EXT(:, 2, 1:13, 1:13)); cond2E(ids, :, :) = []; 
cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); cond1A(ids, :, :) = []; 
cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); cond2A(ids, :,:) = []; 
cond1E = squeeze(out_rsa_EXT(:, 1, :,:)); cond1E(ids, :, :) = []; 
cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); cond2E(ids, :, :) = []; 

diffA = cond1A-cond2A; 
diffE = cond1E-cond2E; 

%diffA = cond1A; 
%diffE = cond1E; 


for subji = 1:size(cond1A, 1)
   diffAB(subji, :) = diag(squeeze(diffA(subji, :, :)));
   diffEB(subji, :) = diag(squeeze(diffE(subji, :, :)));
end       
diffA = diffAB; diffE = diffEB; 



d2pm1	= squeeze(mean(diffA,'omitnan'));
d2pm2	= squeeze(mean(diffE,'omitnan'));
d2pstd1	= std(diffA, 'omitnan');
d2pstd2	= std(diffE, 'omitnan');
se1 = d2pstd1/sqrt(size(diffA, 1));
se2 = d2pstd2/sqrt(size(diffE, 1));


[h p ci ts] = ttest(diffA, diffE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
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
hb = h; hb(h==0) = nan; hb(hb==1) = -.025; 

%times = (-.5:.01:2) + .25;
times = -.25:.1:2.3;
figure(); 
colors2use = brewermap([6],'*Set2')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot(times, hb, LineWidth=8)
set(gca, 'xlim', [-.25 1.7], 'Fontsize', 28);
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 2);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 2);

%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%'ylim', [-.032 .035], 
%title(f2sav, 'Interpreter','none')
exportgraphics(gcf, ['myP.png'], 'Resolution',150)


%% Permutations (2D) 

nPerm = 1000;

nSubj =  size(diffA, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, diffA(:, 6:40), diffE(:, 6:40));
%junts = cat(1, diffA(:, 3:20), diffE(:, 3:20));
%junts = cat(1, cond1(:, 3:22), cond2(:, 3:22));

[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);
    
clear max_clust_sum_perm
for permi = 1:nPerm
    
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);

    cond1P = junts(fakeCondMapping == 0, :);
    cond2P = junts(fakeCondMapping == 1, :);

    diffC = cond1P - cond2P; 
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm









%% plot 2 lines from TG 50ms or 10ms
clear


sub2exc = [19 25 31]; %Always in the 50 subject space


paths = load_paths_EXT; 
f2sav = 'RSA_TMP_C_1-54_1_0_500-10_1_SICSPE-SICSME';
%f2sav = 'RSA_PFC_C_3-54_1_0_500-10_1_SCE-DCE';
load ([ paths.results.rsa f2sav '.mat']);

[out_rsa, ids] = rem_nan_subj_EXT(out_rsa, sub2exc); 

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

% cond1(sub2exc, :) = []; 
% cond2(sub2exc, :) = []; 

d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1));
se2 = d2pstd2/sqrt(size(cond1, 1));

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
else
    tObs = 0; 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
%h = 0; 
hb = h; hb(h==0) = nan; hb(hb==1) = -.005; 

%times = (-.5:.01:2) + .25;
times = -.25:.01:2.25;
figure(); 
%colors2use = brewermap([6],'*Set1')*0.75;
%colors2use = brewermap([6],'*Accent')*0.75;
colors2use = brewermap([6],'Accent')*0.75;

shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(2,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(1,:)}, 1); hold on; 
plot(times, hb, LineWidth=10)
%set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .03], 'Fontsize', 28);
set(gca, 'xlim', [-.25 1.7],'ylim', [-.01 .05], 'Fontsize', 28);
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 3);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 3);

%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%
%title(f2sav, 'Interpreter','none')
%exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)






%% Permutations (2D) 50ms

nPerm = 10000;

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, cond1(:, 51:225), cond2(:, 51:225));
%junts = cat(1, cond1(:, 6:40), cond2(:, 6:40));
%junts = cat(1, cond1(:, 3:22), cond2(:, 3:22));

[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);
    
clear max_clust_sum_perm
for permi = 1:nPerm
    
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);

    cond1P = junts(fakeCondMapping == 0, :);
    cond2P = junts(fakeCondMapping == 1, :);

    diffC = cond1P - cond2P; 
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
%allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm




%% PERMUTATIONS 2D IN LOOP 
clear , clc


nPerm = 10000;



listF2sav = {

'RSA_HPC_C_3-54_1_0_500-100_1_SCE-DCE';
'RSA_AMY_C_3-54_1_0_500-100_1_SCE-DCE';
'RSA_OCC_C_3-54_1_0_500-100_1_SCE-DCE';
'RSA_PFC_C_3-54_1_0_500-100_1_SCE-DCE';
'RSA_OFC_C_3-54_1_0_500-100_1_SCE-DCE';
'RSA_TMP_C_3-54_1_0_500-100_1_SCE-DCE';

'RSA_HPC_C_3-54_1_0_500-100_1_SCA-DCA';
'RSA_AMY_C_3-54_1_0_500-100_1_SCA-DCA';
'RSA_OCC_C_3-54_1_0_500-100_1_SCA-DCA';
'RSA_PFC_C_3-54_1_0_500-100_1_SCA-DCA';
'RSA_OFC_C_3-54_1_0_500-100_1_SCA-DCA';
'RSA_TMP_C_3-54_1_0_500-100_1_SCA-DCA';

'RSA_HPC_C_3-54_1_0_500-100_1_SCT-DCT';
'RSA_AMY_C_3-54_1_0_500-100_1_SCT-DCT';
'RSA_OCC_C_3-54_1_0_500-100_1_SCT-DCT';
'RSA_PFC_C_3-54_1_0_500-100_1_SCT-DCT';
'RSA_OFC_C_3-54_1_0_500-100_1_SCT-DCT';
'RSA_TMP_C_3-54_1_0_500-100_1_SCT-DCT';

'RSA_HPC_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
'RSA_AMY_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
'RSA_OCC_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
'RSA_PFC_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
'RSA_OFC_C_3-54_1_0_500-100_1_SICSPA-SICSMA';
'RSA_TMP_C_3-54_1_0_500-100_1_SICSPA-SICSMA';

'RSA_HPC_C_3-54_1_0_500-100_1_SICSPE-SICSME';
'RSA_AMY_C_3-54_1_0_500-100_1_SICSPE-SICSME';
'RSA_OCC_C_3-54_1_0_500-100_1_SICSPE-SICSME';
'RSA_PFC_C_3-54_1_0_500-100_1_SICSPE-SICSME';
'RSA_OFC_C_3-54_1_0_500-100_1_SICSPE-SICSME';
'RSA_TMP_C_3-54_1_0_500-100_1_SICSPE-SICSME';

'RSA_HPC_C_3-54_1_0_500-100_1_SISVA-DISVA';
'RSA_AMY_C_3-54_1_0_500-100_1_SISVA-DISVA';
'RSA_OCC_C_3-54_1_0_500-100_1_SISVA-DISVA';
'RSA_PFC_C_3-54_1_0_500-100_1_SISVA-DISVA';
'RSA_OFC_C_3-54_1_0_500-100_1_SISVA-DISVA';
'RSA_TMP_C_3-54_1_0_500-100_1_SISVA-DISVA';

'RSA_HPC_C_3-54_1_0_500-100_1_SISVE-DISVE';
'RSA_AMY_C_3-54_1_0_500-100_1_SISVE-DISVE';
'RSA_OCC_C_3-54_1_0_500-100_1_SISVE-DISVE';
'RSA_PFC_C_3-54_1_0_500-100_1_SISVE-DISVE';
'RSA_OFC_C_3-54_1_0_500-100_1_SISVE-DISVE';
'RSA_TMP_C_3-54_1_0_500-100_1_SISVE-DISVE';


};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 4 %:length(listF2sav)
    clearvars -except listF2sav listi t1 paths nPerm allFilesP
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);

    load ([ paths.results.rsa f2sav '.mat']);

    % % % % % % compute real t
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
    se1 = d2pstd1/sqrt(size(cond1, 1));
    se2 = d2pstd2/sqrt(size(cond1, 1));
    
    [h p ci ts] = ttest(cond1, cond2); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        tObs = allSTs(id); 
    else
        tObs = 0; 
    end



    if tObs ~= 0
        nSubj =  size(cond1, 1);
        realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';
        
        junts = cat(1, cond1(:, 3:20), cond2(:, 3:20));
        %junts = cat(1, cond1(:, 3:22), cond2(:, 3:22));
        
        [M,N] = size(realCondMapping);
        rowIndex = repmat((1:M)',[1 N]);
            
        clear max_clust_sum_perm
        for permi = 1:nPerm
            
            [~,randomizedColIndex] = sort(rand(M,N),2);
            newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
            fakeCondMapping = realCondMapping(newLinearIndex);
            fakeCondMapping = fakeCondMapping(:);
        
            cond1P = junts(fakeCondMapping == 0, :);
            cond2P = junts(fakeCondMapping == 1, :);
        
            diffC = cond1P - cond2P; 
            [h p ci ts] = ttest(diffC); 
            h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
            clear allSTs  
            clustinfo = bwconncomp(h);
            for pxi = 1:length(clustinfo.PixelIdxList)
               allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
            end
            
            if exist('allSTs')
                [max2u id] = max(abs(allSTs));
                %[max2u id] = max(allSTs);
                max_clust_sum_perm(permi,:) = allSTs(id); 
            else
                max_clust_sum_perm(permi,:) = 0; 
            end
        
        end
       
    
        %allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
        allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
        allFilesP{listi, 1} = f2sav;
        allFilesP{listi, 2} = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
    
        else
            allFilesP{listi, 1} = f2sav;
            allFilesP{listi, 2} = 1;
        end


end

save ('allFilesP', 'allFilesP')


%% PERMUTATIONS 2D FOR CONTRASTS OF CONTRASTS IN LOOP 
clear, clc
paths = load_paths_EXT; 

nPerm = 1000; 

listF2sav = {



                 {['RSA_TMP_C_3-54_1_0_500-100_1_SICSPA-SICSMA']  ['RSA_TMP_C_3-54_1_0_500-100_1_SICSPE-SICSME'] };
                 {['RSA_OFC_C_3-54_1_0_500-100_1_SICSPA-SICSMA']  ['RSA_OFC_C_3-54_1_0_500-100_1_SICSPE-SICSME'] };
                 {['RSA_HPC_C_3-54_1_0_500-100_1_SICSPA-SICSMA']  ['RSA_HPC_C_3-54_1_0_500-100_1_SICSPE-SICSME'] };
                 {['RSA_AMY_C_3-54_1_0_500-100_1_SICSPA-SICSMA']  ['RSA_AMY_C_3-54_1_0_500-100_1_SICSPE-SICSME'] };
                 {['RSA_PFC_C_3-54_1_0_500-100_1_SICSPA-SICSMA']  ['RSA_PFC_C_3-54_1_0_500-100_1_SICSPE-SICSME'] };

};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths pM nPerm allFilesP
        
    
    f2sav = listF2sav{listi}{1};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_ACQ = out_rsa; 
    f2sav = listF2sav{listi}{2};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_EXT = out_rsa; 

    
    % % % % %  remove hack 
    ids = []; 
    for subji = 1:size(out_rsa, 1)
        %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
        %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
        cond1 = squeeze(out_rsa(subji, 1, :, :)); 
        cond2 = squeeze(out_rsa(subji, 2, :, :)); 
        if cond1(1) == 0
            ids = [ids subji];
        end
    end
    
    % cond1A = squeeze(out_rsa_ACQ(:, 1, 1:13, 1:13)); cond1A(ids, :, :) = []; 
    % cond2A = squeeze(out_rsa_ACQ(:, 2, 1:13, 1:13)); cond2A(ids, :,:) = []; 
    % cond1E = squeeze(out_rsa_EXT(:, 1, 1:13, 1:13)); cond1E(ids, :, :) = []; 
    % cond2E = squeeze(out_rsa_EXT(:, 2, 1:13, 1:13)); cond2E(ids, :, :) = []; 
    cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); cond1A(ids, :, :) = []; 
    cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); cond2A(ids, :,:) = []; 
    cond1E = squeeze(out_rsa_EXT(:, 1, :,:)); cond1E(ids, :, :) = []; 
    cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); cond2E(ids, :, :) = []; 
    
    diffA = cond1A-cond2A; 
    diffE = cond1E-cond2E; 
    
    
    
    [diffA diffE] = rem_half_matrix(diffA, diffE);
    m1 = squeeze(mean(diffA, 'omitnan')); 
    m2 = squeeze(mean(diffE, 'omitnan')); 
    
    [h p ci ts] = ttest(diffA, diffE); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    %h(1:26,1:2) = 0;  % % % no clusters before baseline
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        tObs = allSTs(id); 
    else
        tObs = 0; 
    end



    if tObs ~= 0
        nSubj =  size(cond1, 1);
        realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';
        
        junts = cat(1, cond1(:, 3:20), cond2(:, 3:20));
        %junts = cat(1, cond1(:, 3:22), cond2(:, 3:22));
        
        [M,N] = size(realCondMapping);
        rowIndex = repmat((1:M)',[1 N]);
            
        clear max_clust_sum_perm
        for permi = 1:nPerm
            
            [~,randomizedColIndex] = sort(rand(M,N),2);
            newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
            fakeCondMapping = realCondMapping(newLinearIndex);
            fakeCondMapping = fakeCondMapping(:);
        
            cond1P = junts(fakeCondMapping == 0, :);
            cond2P = junts(fakeCondMapping == 1, :);
        
            diffC = cond1P - cond2P; 
            [h p ci ts] = ttest(diffC); 
            h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
            clear allSTs  
            clustinfo = bwconncomp(h);
            for pxi = 1:length(clustinfo.PixelIdxList)
               allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
            end
            
            if exist('allSTs')
                [max2u id] = max(abs(allSTs));
                %[max2u id] = max(allSTs);
                max_clust_sum_perm(permi,:) = allSTs(id); 
            else
                max_clust_sum_perm(permi,:) = 0; 
            end
        
        end
       
    
        %allAb = max_clust_sum_perm(max_clust_sum_perm > tObs);
        allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
        allFilesP{listi, 1} = f2sav;
        allFilesP{listi, 2} = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
    
        else
            allFilesP{listi, 1} = f2sav;
            allFilesP{listi, 2} = 1;
        end


end

save ('allFilesP', 'allFilesP')

















%% FREQUENCY RESOLVED
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast

clear , clc

listF2sav = {

% 'RSA_AMY_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 'RSA_HPC_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 'RSA_PFC_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 'RSA_OFC_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 'RSA_TMP_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 'RSA_OCC_C_3-54_1_1_500-100_0_SICSPE-SICSME';
% 
% 'RSA_AMY_C_3-54_1_1_500-100_0_SICSPA-SICSMA';
% 'RSA_HPC_C_3-54_1_1_500-100_0_SICSPA-SICSMA';
% 'RSA_PFC_C_3-54_1_1_500-100_0_SICSPA-SICSMA';
% 'RSA_OFC_C_3-54_1_1_500-100_0_SICSPA-SICSMA';
% 'RSA_TMP_C_3-54_1_1_500-100_0_SICSPA-SICSMA';
% 'RSA_OCC_C_3-54_1_1_500-100_0_SICSPA-SICSMA';

% 'RSA_AMY_C_3-54_1_1_500-100_0_SCA-DCA';
% 'RSA_HPC_C_3-54_1_1_500-100_0_SCA-DCA';
% 'RSA_PFC_C_3-54_1_1_500-100_0_SCA-DCA';
% 'RSA_OFC_C_3-54_1_1_500-100_0_SCA-DCA';
% 'RSA_TMP_C_3-54_1_1_500-100_0_SCA-DCA';
% 'RSA_OCC_C_3-54_1_1_500-100_0_SCA-DCA';
% 
%'RSA_AMY_C_3-54_1_1_500-100_0_SCE-DCE';
% 'RSA_HPC_C_3-54_1_1_500-100_0_SCE-DCE';
% 'RSA_PFC_C_3-54_1_1_500-100_0_SCE-DCE';
% 'RSA_OFC_C_3-54_1_1_500-100_0_SCE-DCE';
% 'RSA_TMP_C_3-54_1_1_500-100_0_SCE-DCE';
% 'RSA_OCC_C_3-54_1_1_500-100_0_SCE-DCE';
% 
% 'RSA_AMY_C_3-54_1_1_500-100_0_SCT-DCT';
% 'RSA_HPC_C_3-54_1_1_500-100_0_SCT-DCT';
% 'RSA_PFC_C_3-54_1_1_500-100_0_SCT-DCT';
% 'RSA_OFC_C_3-54_1_1_500-100_0_SCT-DCT';
% 'RSA_TMP_C_3-54_1_1_500-100_0_SCT-DCT';
% 'RSA_OCC_C_3-54_1_1_500-100_0_SCT-DCT';
% 
% 
% 'RSA_AMY_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 'RSA_HPC_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 'RSA_PFC_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 'RSA_OFC_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 'RSA_TMP_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 'RSA_OCC_C_3-54_0_1_500-100_0_SICSPE-SICSME';
% 
% 'RSA_AMY_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
% 'RSA_HPC_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
% 'RSA_PFC_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
% 'RSA_OFC_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
% 'RSA_TMP_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
% 'RSA_OCC_C_3-54_0_1_500-100_0_SICSPA-SICSMA';

'RSA_AMY_C_3-54_0_1_500-100_0_SCA-DCA';
'RSA_HPC_C_3-54_0_1_500-100_0_SCA-DCA';
'RSA_PFC_C_3-54_0_1_500-100_0_SCA-DCA';
'RSA_OFC_C_3-54_0_1_500-100_0_SCA-DCA';
'RSA_TMP_C_3-54_0_1_500-100_0_SCA-DCA';
'RSA_OCC_C_3-54_0_1_500-100_0_SCA-DCA';

'RSA_AMY_C_3-54_0_1_500-100_0_SCE-DCE';
'RSA_HPC_C_3-54_0_1_500-100_0_SCE-DCE';
'RSA_PFC_C_3-54_0_1_500-100_0_SCE-DCE';
'RSA_OFC_C_3-54_0_1_500-100_0_SCE-DCE';
'RSA_TMP_C_3-54_0_1_500-100_0_SCE-DCE';
'RSA_OCC_C_3-54_0_1_500-100_0_SCE-DCE';



};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
        
    f2sav       = listF2sav{listi}; 
    cfg         = getParams_EXT(f2sav);
    freqs2comp  = cfg.freqs; 
    
    ALLPOW = load ([paths.results.POWfromRT cfg.powF2load]); 
    
    
    for subji = 1:length(ALLPOW.POW)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        cfg.oneListPow = ALLPOW.POW{subji, 1};
        cfg.oneListIds = ALLPOW.POW{subji, 2};
        
        if ~isempty(cfg.oneListPow)
            out_contrasts = create_contrasts_EXT(cfg);
            ids{subji,:} = out_contrasts.allIDs;
    
            for freqi = 1:length(freqs2comp)
                cfg.freqs = freqs2comp(freqi);
                out_rsa(subji, freqi, :, :) = rsa_EXT(out_contrasts, cfg);
            end
        end
        
    end

    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%% plot TF
clear, clc
paths = load_paths_EXT; 
 
f2sav =  'RSA_TMP_C_3-54_0_1_500-100_0_SICSPA-SICSMA';
sub2exc = [];


load ([ paths.results.rsa f2sav '.mat']);

out_rsa(sub2exc, :, :,:) = []; 
ids = rem_nan_subj_EXT(out_rsa); 

cond1 = squeeze(out_rsa(:, :, 1, 1:23)); 
cond2 = squeeze(out_rsa(:, :, 2, 1:23)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
diff = cond1-cond2; 

%[cond1 cond2] = rem_half_matrix(cond1, cond2);

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
%h(1:251,1:26) = 0;  % % % no clusters before baseline

clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    %[max2u id] = max((allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


clim = [-.03 .03];
plot_TG_map(m1, m2, h, t, f2sav, clim)
exportgraphics(gcf, ['myP.png'], 'Resolution',150)












































%%
