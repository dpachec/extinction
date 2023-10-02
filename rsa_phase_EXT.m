%%
%% Temporal RSA
clear

paths = load_paths_EXT; 
file2load = ['TR_' 'Amygdala' '_V']; 
%file2load = ['allS_' 'Hippocampus' '_C']; 
%file2load = ['allS_' 'orbitofrontal' '_C']; 
%file2load = ['allS_' 'superiorfrontal' '_C']; 

load ([paths.results.traces file2load]); 


%% 

%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_trials/noTrials_win-width_mf
clc
clearvars -except ALLEEG paths file2load

f2sav =  'T_nan_1_0_0_50-10_1_SCA-DCA';
%f2sav = 'T_nan_0_0_0_50-10_1_DISVE-DIDVE';

cfg = getParams_EXT(f2sav);

t1 = datetime; 
for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        
        cfg.oneListIds = Ev2; 
        
        EEG = add_EEGLAB_fields(EEG); 
        %EEG = pop_eegfiltnew (EEG, .1,30); %low pass filter up to 30Hz (kunz 2019)
%         d2p = squeeze(EEG.data(1, 1:1000, 1))
%         figure()
%         plot(d2p)

        EEG = normalize_baseline_EXT(EEG, [2501:3000]); 
        EEG = downsample_EEG_EXT(EEG); 
        cfg.oneListTraces = permute(EEG.data(:, 251:500,:), [3 1 2]); 
        %cfg.oneListTraces = permute(EEG.data(:, 2501:5000,:), [3 1 2]); 
        cfg.tyRSA = 'tRSA'; 
        out_contrasts = create_contrasts_EXT(cfg);
        
        out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
        
    end

end


mkdir ([paths.results.rsa]);
ids = out_contrasts.allIDs; 
save([ paths.results.rsa f2sav '_' file2load '.mat'], 'out_rsa', 'ids');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%%
clear
paths = load_paths_EXT; 
file2load = ['TR_' 'Amygdala' '_C']; 
f2sav = 'T_nan_0_0_0_50-10_1_DISVE-DIDVE';
load ([ paths.results.rsa f2sav '_' file2load '.mat']);



%% remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)

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


% % % remove half of the matrix
for subji = 1:size(cond1, 1)
    rdm2Tril = squeeze(cond1(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond1TR(subji, :, :) = rdm2Tril;

    rdm2Tril = squeeze(cond2(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond2TR(subji, :, :) = rdm2Tril;
end


m1 = squeeze(mean(cond1TR, 'omitnan')); 
m2 = squeeze(mean(cond2TR, 'omitnan')); 

[h p ci ts] = ttest(cond1TR, cond2TR); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

[max2u id] = max(abs(allSTs));
tObs = allSTs(id); 


h = zeros(size(cond1TR, 2),size(cond1TR, 2)); 
h(clustinfo.PixelIdxList{id}) = 1;


figure(); tiledlayout(1,3);
nexttile
%imagesc(m1);  axis square
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; colorbar
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-.04 .04])

nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square;hold on; colorbar 
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-.04 .04])

nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-3 3])


axesHandles = findall(0, 'type', 'axes');
%set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 150], 'ylim', [1 150]); 
set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', []); 
%colorbar
exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)


%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, cond1TR(:, 6:15, 6:15), cond2TR(:, 6:15, 6:15));

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
clear, close all
paths = load_paths_EXT; 

c2u = 'V';
sROI = {'Amygdala'}; 

%sROI = {'superiorfrontal'}; 

%sROI = {'superiorfrontal' 'rostralmiddlefrontal' 'anteriorcingulate' 'posteriorcingulate' 'precentral' 'caudalmiddlefrontal'}; % case sensitive 

 %sROI = { 'inferiortemporal' 'middletemporal' 'superiortemporal' 'bankssts' 'fusiform' 'temporalpole' ...
 %             'lateraloccipital' 'lingual' 'parahippocampal' 'cuneus' 'pericalcarine' };

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';



for subji = 1:length(allsubs)
    subji
    sub = allsubs{subji}; 
    cd([ paths.iEEG])
    load ([sub '_iEEG.mat']);

    EEG = remove_elec_EXT_manually(EEG, subji); %thres channels is 1/5 of 192 = 38
    chansLab = {EEG.chanlocs.fsLabel}';
    selChans = contains(chansLab, sROI);

    
    if find(selChans)
        
        EEG.chanlocs = EEG.chanlocs(selChans);
        EEG.data = EEG.data(selChans, :); %contains nans
        
        %epoch data
        Ev = [{EEG.event.type}]'; 
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        clen = cellfun(@length, Ev1); 
        EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
        Ev2 = cat(1, Ev1{:});
       
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %paris subjects have a space in the last character of the event WHY??
        ids = strcmp(Ev2(:, 10), c2u); 
        
        EEG.event = EEG.event(ids);
        
        Ev2 = Ev2(ids,:);  %store ids
        allEv{subji,:} = Ev2; 
        
        % extract phase time-series
        EEG = add_EEGLAB_fields(EEG);
        EEG = remove_elec_EXT(EEG, 50); %thres channels is 1/5 of 192 = 38

        nChans = size(EEG.data, 1); 
        %phaTS = zeros(192, nChans, 30, 700); 
        clear phaTS; 
    
        for freqi = 3:30
    
            EEG1 = pop_eegfiltnew (EEG, freqi,freqi);
            %EEG.phaTS = angle(hilbert(EEG.data(2,:)));
            %eegplot(EEG.phaTS, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1, 'events', EEG.event);
    
            %epoch data and markers
            [EEG1 id2check] = pop_epoch( EEG1, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
            EEG1 = normalize_baseline_EXT(EEG1, [2501:3000]);  
           
            nChans = size(EEG1.data, 1); nTrials = size(EEG1.data, 3); 
            for chani = 1:nChans
                for triali = 1:nTrials
                    data        = squeeze(EEG1.data(chani, 2501:4700, triali)); % until 1s
                    dataHA      = angle(hilbert(data));
                    dataHADS      = downsample(dataHA, 10);
                    phaTS(triali, chani, freqi,:) = dataHADS; 
                end
            end
        end
        allPHA{subji,:} = phaTS; 
    end
end

sROI = char(join(sROI, '_'));
mkdir(paths.results.phases)
filename = [paths.results.phases 'allPHA_3-30Hz_'  char(sROI) '_' c2u];
save(filename, 'allPHA', 'allEv', '-v7.3');

cd (paths.github)
% 
%% Load Phase time-series
clear

paths = load_paths_EXT; 
file2load = ['allPHA_3-30Hz_Amygdala_V']; 
%file2load = ['allS_' 'Hippocampus' '_C']; 
%file2load = ['allS_' 'orbitofrontal' '_C']; 
%file2load = ['allS_' 'superiorfrontal' '_C']; 

load ([paths.results.phases file2load]); 

%% quantify nans
for subji = 1
    phaTS = allPHA{subji};
    
    for chani = 1%:size(phaTS, 2)
        for triali = 1:20 %size(phaTS, 1)
            d2p = squeeze(phaTS(triali, chani, :, :)); 
            %figure()
            %imagesc(d2p)
            figure()
            plot(d2p(8,:))
        end
    end
end


%% remove 1Hz and 2Hz
for subji = 1:47
    phaTS = allPHA{subji};
    if ~isempty(phaTS)
        for chani = 1:size(phaTS, 2)
            for triali = 1:size(phaTS, 1)
                d2p = squeeze(phaTS(triali, chani, :, :)); 
                d2p = d2p(3:30, :);
                phaTS2(triali, chani, :, :) = d2p; 
            end
        end
        allPHA{subji} = phaTS2; 
    else
        allPHA{subji} = []; 
    end
end


%% Phases as patterns based RSA
%freqs_avTimeFeatVect_freqResolv(0-1)_trials/noTrials_win-width_mf
clc
clearvars -except allPHA allEv paths

f2sav = 'PHA_1-28_0_1_0_50-10_1_SCA-DCA'; 

cfg = getParams_EXT(f2sav);

t1 = datetime; 
for subji = 1:length(allPHA)
    
    phaTS = allPHA{subji};
   
    
    if ~isempty(phaTS)
        Ev2   = allEv{subji}; 
        ids = strcmp(Ev2(:, 10), 'V'); 
        Ev2 = Ev2(ids,:); 

        cfg.tyRSA = 'pRSA'; 
        cfg.oneListIds = Ev2; 
        cfg.oneListPow = phaTS(:, :, : ,:); 
        
        out_contrasts = create_contrasts_EXT(cfg);
        
        tic
        %out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
        %out_rsa(subji, :, :, :) = rsa_EXT2(out_contrasts, cfg);
        out_rsa(subji, :, :, :) = rsa_EXT3(out_contrasts, cfg);
        %out_rsa(subji, :, :, :) = rsa_EXT5(out_contrasts, cfg);
        toc
    
    end

end


mkdir ([paths.results.phases]);
save([ paths.results.phases f2sav '_' 'Amygdala_V.mat'], 'out_rsa', 'Ev2');

t2 = datetime; 
etime(datevec(t2), datevec(t1))

%% load out_rsa
clearvars 
f2sav = 'PHA_1-28_0_1_0_50-10_1_SCA-DCA'; 
paths = load_paths_EXT; 
load ([ paths.results.phases f2sav '_Amygdala_V.mat']);


%% remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)

    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(100) == 0 | isnan(cond1(100))
        ids = [ids subji];
    end

end


%%
cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 


% % % remove half of the matrix
for subji = 1:size(cond1, 1)
    rdm2Tril = squeeze(cond1(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond1TR(subji, :, :) = rdm2Tril;

    rdm2Tril = squeeze(cond2(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond2TR(subji, :, :) = rdm2Tril;
end


m1 = squeeze(mean(cond1TR, 'omitnan')); 
m2 = squeeze(mean(cond2TR, 'omitnan')); 

[h p ci ts] = ttest(cond1TR, cond2TR); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

[max2u id] = max(abs(allSTs));
tObs = allSTs(id); 


%h = zeros(size(cond1TR, 2),size(cond1TR, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


figure(); tiledlayout(1,3);
nexttile
%imagesc(m1);  axis square
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; colorbar
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-.03 .03])

nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square;hold on; colorbar 
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-.03 .03])

nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
plot(get(gca,'xlim'), [5 5],'k', 'linewidth', 1); plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-4 4])


axesHandles = findall(0, 'type', 'axes');
%set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 150], 'ylim', [1 150]); 
set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', []); 
%colorbar
exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)




%% Plot time-frequency 


cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

[max2u id] = max(abs(allSTs));
tObs = allSTs(id); 


%h = zeros(size(cond1TR, 2),size(cond1TR, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


figure(); tiledlayout(1,3);
nexttile
%imagesc(m1);  axis square
contourf( m1, 50, 'linecolor', 'none'); axis square; hold on; colorbar
plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
%set(gca, 'clim', [-.03 .03])

nexttile
contourf( m2, 50, 'linecolor', 'none'); axis square;hold on; colorbar 
plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
%set(gca, 'clim', [-.03 .03])

nexttile
contourf( t, 50, 'linecolor', 'none'); axis square; hold on; colorbar
contour( h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
plot([5 5], get(gca,'ylim'),'k', 'linewidth', 1); 
%plot(get(gca,'xlim'), [25 25],'k', 'linewidth', 1); plot([25 25], get(gca,'ylim'),'k', 'linewidth', 1); 
set(gca, 'clim', [-4 4])


axesHandles = findall(0, 'type', 'axes');
%set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', [], 'xlim', [1 150], 'ylim', [1 150]); 
set(axesHandles,'xtick', [], 'xticklabel', [], 'ytick', [], 'yticklabel', []); 
%colorbar
exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)








%% PLV inter area (PFV-VVS)
%% particular time period and band

clearvars 


f2u = [3 8];
tP = 2401:3200;
%tP = 1201:2000;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)


    t_vvs = squeeze(c_vvs.oneListTraces(:,tP,:));
    EEG = []; EEG.data = t_vvs; EEG.trials = size(t_vvs, 3); EEG.srate=1000; EEG.nbchan=size(t_vvs, 1); EEG.pnts=size(t_vvs,2);EEG.event=[];
    EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
    data_vvs        = squeeze(EEG.data); 
    dataHA_vvs      = angle(hilbert(data_vvs));
    t_pfc = squeeze(c_pfc.oneListTraces(:,tP,:));
    EEG = []; EEG.data = t_pfc; EEG.trials = size(t_pfc, 3); EEG.srate=1000; EEG.nbchan=size(t_pfc, 1); EEG.pnts=size(t_pfc,2);EEG.event=[];
    EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
    data_pfc        = squeeze(EEG.data); 
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        parfor chanj = 1:size(c_pfc.chanNames,1)
            diffPha = angle(hilbert(squeeze(data_vvs(chani, :,:)))) - angle(hilbert(squeeze(data_pfc(chanj, :,:))));
            PLV2U(:, chani, chanj) = abs(mean(exp(1i*(diffPha))));
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


save([paths.results.PLV 'time_PLV_ALL_' num2str(f2u(1)) '-' num2str(f2u(2)) '_' num2str(tP(1)) '-' num2str(tP(end))], 'PLV_ALL');




%% process PLV_ALL (across time) 

clear 
paths = load_paths_WM('vvs')
load ([paths.results.PLV 'time_PLV_ALL_3-8_1201-2000'])
PLV_ALL_Baseline = PLV_ALL;
load ([paths.results.PLV 'time_PLV_ALL_3-8_2401-3200'])

clear SI_TR MI_TR

for subji = 1:10

    allPLV = PLV_ALL{subji, 1};
     allPLVB = PLV_ALL_Baseline{subji, 1}; 
     mT = mean(allPLVB);
     stdT = std(allPLVB);
     allPLV = bsxfun(@rdivide, allPLV - mT, stdT); 

    allIDs = PLV_ALL{subji, 2};
    ids0 = cellfun(@(x) strsplit(x), allIDs, 'un', 0);
    ids1 = cell2mat(cellfun(@(x) double(string(x(1))), ids0, 'un', 0));
    ids2 = cell2mat(cellfun(@(x) double(string(x(2))), ids0, 'un', 0));

    
    SI_TR{subji,1} = allPLV(ids1 == 7 & ids2 ~= 4,:,:,:); 
    MI_TR{subji,1} = allPLV(ids1 == 7 & ids2 == 4,:,:,:); 

end

%% plot 

d2pSI = cellfun(@(x) squeeze(mean(mean(mean(x)))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(mean(x)))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 
%c1 = logit(c1)
%c2 = logit(c2)
md2pSI = mean(c1)
md2pMI = mean(c2)

[h p ci ts] = ttest(c1, c2)

%% 2Bar 
data.data = [c1 c2]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
%set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.32 .45] );
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.05 .25] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


%% collect all with baseline and plot together

clear 

b2u = {'3-8' '9-12' '13-29' '30-75' '75-150'};

for bandi = 1:5

    paths = load_paths_WM('vvs');
    load ([paths.results.PLV 'time_PLV_ALL_' b2u{bandi} '_1201-2000'])
    PLV_ALL_Baseline = PLV_ALL;
    load ([paths.results.PLV 'time_PLV_ALL_' b2u{bandi} '_2001-2800'])
    
    clear SI_TR MI_TR
    
    for subji = 1:10
    
        allPLV = PLV_ALL{subji, 1};
        allPLVB = PLV_ALL_Baseline{subji, 1}; 
        mT = mean(allPLVB);
        stdT = std(allPLVB);
        allPLV = bsxfun(@rdivide, allPLV - mT, stdT); 
    
        allIDs = PLV_ALL{subji, 2};
        ids0 = cellfun(@(x) strsplit(x), allIDs, 'un', 0);
        ids1 = cell2mat(cellfun(@(x) double(string(x(1))), ids0, 'un', 0));
        ids2 = cell2mat(cellfun(@(x) double(string(x(2))), ids0, 'un', 0));
    
        
        SI_TR{subji,1} = allPLV(ids1 == 7 & ids2 ~= 4,:,:,:); 
        MI_TR{subji,1} = allPLV(ids1 == 7 & ids2 == 4,:,:,:); 
    

    end

    d2pSI = cellfun(@(x) squeeze(mean(mean(mean(x)))), SI_TR, 'un', 0);
    d2pMI = cellfun(@(x) squeeze(mean(mean(mean(x)))), MI_TR, 'un', 0);
    c1(:, bandi) = cell2mat(d2pSI')'; 
    c2(:, bandi) = cell2mat(d2pMI')'; 


end



%% 
mC1 = squeeze(mean(c1));
stdC1 = std(c1, [], 1); 
seC1 = stdC1 / sqrt(10);

mC2 = squeeze(mean(c2));
stdC2 = std(c2, [], 1); 
seC2 = stdC2 / sqrt(10);

shadedErrorBar(1:5, mC1, seC1, 'r', 1); hold on; 
shadedErrorBar(1:5, mC2, seC2, 'b', 1); hold on; 

[h p ci ts] = ttest(c1, c2)

%% 6Bar 
data.data = [c1(:,1) c2(:,1) c1(:,2) c2(:,2) c1(:,3) c2(:,3) c1(:,4) c2(:,4) c1(:,5) c2(:,5)]; 
%data.data = [c1  c2]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data,1);
hb = plot ([1:10], data.data', 'linestyle','none'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h, 'Color','k','linestyle','none', 'lineWidth', 2);
%set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.32 .45] );
set(gca,'XTick',[1:10 ],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 11], 'ylim', [-.2 .35] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);







%% CONNECTIVITY Across trials
% % % % PLV PLI and COHERENCE

clearvars 

f2u = [3 8];

paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)
    
    clear PLV_SI PLV_MI PLI_SI PLI_MI PLV_SI_M2 PLV_MI_M2 PLI_SI_M2 PLI_MI_M2 ...
    WPLI_SI WPLI_MI COH_SI COH_MI ANG_SI ANG_MI dPHA_SI dPHA_MI

    for chani = 1:size(c_vvs.chanNames,1)
        parfor chanj = 1:size(c_pfc.chanNames,1)
            id2u = cellfun(@(x) strsplit(x, ' '), c_vvs.oneListIds_c, 'un', 0);
            id2u = cat(1, id2u{:});
            id2u_SI = strcmp(id2u(:, 1), '7') & ~strcmp(id2u(:, 2), '4'); 
            id2u_MI = strcmp(id2u(:, 1), '7') & strcmp(id2u(:, 2), '4'); 

            % % 
            t_vvsSI = squeeze(c_vvs.oneListTraces(chani,:,id2u_SI));
            t_pfcSI = squeeze(c_pfc.oneListTraces(chanj,:,id2u_SI));
            t_vvsMI = squeeze(c_vvs.oneListTraces(chani,:,id2u_MI));
            t_pfcMI = squeeze(c_pfc.oneListTraces(chanj,:,id2u_MI));            
            
            %match trial numbers
            % get condition with more trials 
            nTrC1 = size(t_vvsSI,2); nTrC2 = size(t_vvsMI, 2); 
            id2uSample = randperm(nTrC2);
            if nTrC1 > nTrC2
                t_vvsSI = t_vvsSI(:, id2uSample);
                t_pfcSI = t_pfcSI(:, id2uSample);
            elseif nTrC2 > nTrC1
                t_vvsMI = t_vvsMI(:, id2uSample);
                t_pfcMI = t_pfcMI(:, id2uSample);
            end

            % % % SINGLE ITEM TRIALS
            EEG = []; 
            EEG.data    = t_vvsSI;
            EEG.trials  = size(t_vvsSI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_vvsSI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_vvs        = hilbert(squeeze(EEG.data)); 
            dataHA_vvs      = angle(data_vvs);
            EEG = []; 
            EEG.data    = t_pfcSI;
            EEG.trials  = size(t_pfcSI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_pfcSI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_pfc        = hilbert(squeeze(EEG.data)); 
            dataHA_pfc      = angle(data_pfc);
            diffPha = dataHA_vvs - dataHA_pfc;
            PLV_SI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI_SI(chani, chanj, :) = abs(mean(sign(imag(exp(1i*diffPha))),2)); %PLI
            cdd = data_vvs .* conj(data_pfc);% cross-spectral density
            cdi = imag(cdd);
            PLV_SI_M2(chani, chanj, :) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI_SI_M2(chani, chanj, :) = abs(mean(sign(imag(cdd)),2));
            WPLI_SI(chani, chanj, :) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
            % % compute coherence
            spec1 = mean(data_vvs.*conj(data_vvs),2);
            spec2 = mean(data_pfc.*conj(data_pfc),2);
            specX = abs(mean(data_vvs.*conj(data_pfc),2)).^2;
            COH_SI(chani, chanj, :) = specX./ (spec1.*spec2);
            ANG_SI(chani, chanj, :) = angle(mean(exp(1i*(diffPha)), 2)); 
            dPHA_SI(chani, chanj, :,:) = diffPha; 
            dPHA2_SI(chani, chanj, :,:) = diffPha2; 


            % % % % MULTI ITEM TRIALS
            EEG = []; 
            EEG.data    = t_vvsMI;
            EEG.trials  = size(t_vvsMI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_vvsMI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_vvs        = hilbert(squeeze(EEG.data)); 
            dataHA_vvs      = angle(data_vvs);
            EEG = []; 
            EEG.data    = t_pfcMI;
            EEG.trials  = size(t_pfcMI, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_pfcMI,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_pfc        = hilbert(squeeze(EEG.data)); 
            dataHA_pfc      = angle(data_pfc);
            diffPha = dataHA_vvs - dataHA_pfc;
            PLV_MI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI_MI(chani, chanj, :) = abs(mean(sign(imag(exp(1i*diffPha))),2)); %PLI
            cdd = data_vvs .* conj(data_pfc);% cross-spectral density
            cdi = imag(cdd);
            PLV_MI_M2(chani, chanj, :) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI_MI_M2(chani, chanj, :)  = abs(mean(sign(imag(cdd)),2));            
            WPLI_MI(chani, chanj, :) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
            % % compute coherence
            spec1 = mean(data_vvs.*conj(data_vvs),2);
            spec2 = mean(data_pfc.*conj(data_pfc),2);
            specX = abs(mean(data_vvs.*conj(data_pfc),2)).^2;
            COH_MI(chani, chanj, :) = specX./ (spec1.*spec2);
            ANG_MI(chani, chanj, :) = angle(mean(exp(1i*(diffPha)), 2)); 
            dPHA_MI(chani, chanj, :,:) = diffPha; 

        end
    end
    CON_ALL{subji,1} = PLV_SI; %pfc or vvs 
    CON_ALL{subji,2} = PLV_MI; %pfc or vvs 
    CON_ALL{subji,3} = PLI_SI; %pfc or vvs 
    CON_ALL{subji,4} = PLI_MI; %pfc or vvs 
    CON_ALL{subji,5} = PLV_SI_M2; %pfc or vvs 
    CON_ALL{subji,6} = PLV_MI_M2; %pfc or vvs 
    CON_ALL{subji,7} = PLI_SI_M2; %pfc or vvs 
    CON_ALL{subji,8} = PLI_MI_M2; %pfc or vvs 
    CON_ALL{subji,9} = WPLI_SI; %pfc or vvs 
    CON_ALL{subji,10} = WPLI_MI; %pfc or vvs 
    CON_ALL{subji,11} = COH_SI; %pfc or vvs 
    CON_ALL{subji,12} = COH_MI; %pfc or vvs 
    CON_ALL{subji,13} = ANG_SI; %pfc or vvs 
    CON_ALL{subji,14} = ANG_MI; %pfc or vvs 
    CON_ALL{subji,15} = dPHA_SI; %pfc or vvs 
    CON_ALL{subji,16} = dPHA_MI; %pfc or vvs 

    

end


save([paths.results.PLV 'trials_CON_ALL_' num2str(f2u(1)) '-' num2str(f2u(2)) 'Hz'], 'CON_ALL');



toc






%% Process across trials 

clear 
paths = load_paths_WM('vvs')
load ([paths.results.PLV 'trials_CON_ALL_3-8Hz'])



%% 


SI_TR = CON_ALL(:, 1);
MI_TR = CON_ALL(:, 2);

d2pSI = cellfun(@(x) squeeze(mean(mean(x))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(x))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 


% % % normalize to the baseline
% mC1 = mean(c1(:, 1000:2000), 2);
% stdC1 = std(c1(:, 1000:2000), [], 2);
% c1 = bsxfun(@rdivide, c1 - mC1, stdC1); 
% mC2 = mean(c2(:, 1000:2000), 2);
% stdC2 = std(c2(:, 1000:2000), [], 2);
% c2 = bsxfun(@rdivide, c2 - mC2, stdC2); 


md2pSI = mean(c1);
md2pMI = mean(c2);

[h p ci ts] = ttest(c1, c2);


%% Plot
d2SI = cell2mat(d2pSI')';
d2MI = cell2mat(d2pMI')';

mSI = squeeze(mean(d2SI));
mMI = squeeze(mean(d2MI));

times = -2:0.001:6.9999;
hb = h; hb(hb==0) = nan; hb(hb==1) = 0.1; 
figure;
plot(times, mSI, 'r', 'LineWidth', 2); hold on
plot(times, mMI, 'b', 'LineWidth', 2)
plot(times, hb, 'k', 'LineWidth', 8)
legend({'SI' 'MI'})

set(gca, 'FontSize', 14, 'xlim', [-1 3.5], 'ylim', [.1 .2])

%% check for the time period of network fit

c1R = mean(c1(:, 2400:3200), 2)
c2R = mean(c2(:, 2400:3200), 2)
[h p ci ts] = ttest(c1R, c2R);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);
%[pval, F] = circ_htest(c1R, c2R)



%% 2Bar 
data.data = [c1R c2R]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.27 .64] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);







%% COMPUTE AVERAGE PHASE LAG
clearvars -except CON_ALL

for subji = 1:10
    diffPha_SI = CON_ALL{subji, 15}; 
    diffPha_MI = CON_ALL{subji, 16}; 
    countA = 1; 
    clear allSI allMI
    for chani = 1:size(diffPha_SI, 1)
        for chanj = 1:size(diffPha_SI, 2)

            diffPHASIH = squeeze(diffPha_SI(chani, chanj, :, :)); 
            diffPhaSI_2 = angle(mean(exp(1i*(diffPHASIH)), 2)); 
            diffPhaSI_3{subji,:}(chani, chanj, :) = angle(mean(exp(1i*(diffPhaSI_2)))); 
            allSI(countA, :) = angle(mean(exp(1i*(diffPhaSI_2)))); 


            diffPHAMIH = squeeze(diffPha_MI(chani, chanj, :, :)); 
            diffPhaMI_2 = angle(mean(exp(1i*(diffPHAMIH)), 2)); 
            diffPhaMI_3{subji,:}(chani, chanj, :) = angle(mean(exp(1i*(diffPhaMI_2)))); 
            allMI(countA, :) = angle(mean(exp(1i*(diffPhaMI_2)))); 

            countA = countA+1; 
            
           
        end
    end

        figure()
        tiledlayout(1, 2)
        nexttile
        histogram(allSI, 10); 
        set(gca, 'xlim', [-4 4])
        nexttile
        histogram(allMI, 10); 
        set(gca, 'xlim', [-4 4])
        exportgraphics(gcf, [num2str(subji), '_.png'], 'Resolution', 150)
end
    
%% 
figure; 
histogram(allMI)


%% 
for subji = 1:10

%     allSI = diffPhaSI_3{subji}; 
%     figure()
%     histogram(allSI)


    allMI = diffPhaMI_3{subji}; 
    figure()
    histogram(allMI)

end





%%

SI_TR = diffPhaSI_3; 
MI_TR = diffPhaMI_3; 
d2pSI = cellfun(@(x) squeeze(mean(mean(x))), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(x))), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI')'; 
c2 = cell2mat(d2pMI')'; 


%% 2Bar 
data.data = [c1 ; c2]'; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
%set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.32 .45] );
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.25 .25] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


%% quantify in how many channels angle differences are different from zero 


for subji = 1:10
    
    dPHA_SI = squeeze(mean(CON_ALL{subji, 15}(:, :, 2400:3200,:), 3)); 
    dPHA_MI = squeeze(mean(CON_ALL{subji, 16}(:, :, 2400:3200,:), 3)); 

    for chani = 1:size(dPHA_SI, 1)

        for chanj = 1:size(dPHA_SI, 2)

            dPHA_tr_SI = dPHA_SI(chani, chanj, :); 
            allHs_SI{subji}(chani, chanj,:) = ttest(dPHA_tr_SI);

            dPHA_tr_MI = dPHA_MI(chani, chanj, :); 
            allHs_MI{subji}(chani, chanj,:) = ttest(dPHA_tr_MI);

        end

    end

    figure();
    tiledlayout(1, 2)
    nexttile
    imagesc(allHs_SI{subji});axis square
    title('Single-item')
    nexttile
    imagesc(allHs_MI{subji});axis square
    title('Multi-item')

    filename = num2str(subji)
    exportgraphics(gcf, [filename '.png'], 'Resolution', 150)
end







%% CONNECTIVITY Across trials for all trials (no different conditions)


clearvars 

f2u = [3 8];

paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)
    
    clear PLV PLI PLV_M2 PLI_M2 WPLI COH ANG dPHA

    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            id2u = cellfun(@(x) strsplit(x, ' '), c_vvs.oneListIds_c, 'un', 0);
            id2u = cat(1, id2u{:});
            id2u = strcmp(id2u(:, 1), '7') & ~strcmp(id2u(:, 2), '4'); 

            % % 
            t_vvs = squeeze(c_vvs.oneListTraces(chani,:,id2u));
            t_pfc = squeeze(c_pfc.oneListTraces(chanj,:,id2u));
                        
            
            EEG = []; 
            EEG.data    = t_vvs;
            EEG.trials  = size(t_vvs, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_vvs,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_vvs        = hilbert(squeeze(EEG.data)); 
            dataHA_vvs      = angle(data_vvs);
            EEG = []; 
            EEG.data    = t_pfc;
            EEG.trials  = size(t_pfc, 2); EEG.srate   = 1000; EEG.nbchan  = 1; 
            EEG.pnts = size(t_pfc,1);EEG.event   = [];
            EEG         = pop_eegfiltnew (EEG, f2u(1),f2u(2));
            data_pfc        = hilbert(squeeze(EEG.data)); 
            dataHA_pfc      = angle(data_pfc);
            diffPha = dataHA_vvs - dataHA_pfc;
            PLV(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI(chani, chanj, :) = abs(mean(sign(imag(exp(1i*diffPha))),2)); %PLI
            cdd = data_vvs .* conj(data_pfc);% cross-spectral density
            cdi = imag(cdd);
            PLV_M2(chani, chanj, :) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(chani, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
            PLI_M2(chani, chanj, :) = abs(mean(sign(imag(cdd)),2));
            WPLI(chani, chanj, :) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
            % % compute coherence
            spec1 = mean(data_vvs.*conj(data_vvs),2);
            spec2 = mean(data_pfc.*conj(data_pfc),2);
            specX = abs(mean(data_vvs.*conj(data_pfc),2)).^2;
            COH(chani, chanj, :) = specX./ (spec1.*spec2);
            ANG(chani, chanj, :) = angle(mean(exp(1i*(diffPha)), 2)); %PLV
            dPHA(chani, chanj, :,:) = diffPha; 


        end
    end
    CON_ALL{subji,1} = PLV; %pfc or vvs 
    CON_ALL{subji,2} = PLI; %pfc or vvs 
    CON_ALL{subji,3} = PLV_M2; %pfc or vvs 
    CON_ALL{subji,4} = PLI_M2; %pfc or vvs 
    CON_ALL{subji,5} = WPLI; %pfc or vvs 
    CON_ALL{subji,6} = COH; %pfc or vvs 
    CON_ALL{subji,7} = ANG; %pfc or vvs 
    CON_ALL{subji,8} = dPHA; %pfc or vvs 
    CON_ALL{subji,9} = c_vvs.oneListIds_c(id2u); %pfc or vvs 

    

end


save([paths.results.PLV 'allT_CON_ALL_SI' num2str(f2u(1)) '-' num2str(f2u(2)) 'Hz'], 'CON_ALL');



toc



%% Process across trials 

clear 
paths = load_paths_WM('vvs')
load ([paths.results.PLV 'allT_CON_ALL_SI3-8Hz'])


%% quantify in how many channels angle differences are different from zero 


for subji = 1:10 
    
    dPHA_SI = squeeze(mean(CON_ALL{subji, 8}(:, :, 2400:3200,:), 3)); 

    for chani = 1:size(dPHA_SI, 1)

        for chanj = 1:size(dPHA_SI, 2)

            dPHA_tr = dPHA_SI(chani, chanj, :); 
            allHs{subji}(chani, chanj,:) = ttest(dPHA_tr);



        end

    end


end

%% extract single trial metric of phase differences


clear disT2M
for subji = 1:10

    phaDiff = CON_ALL{subji, 8};
    
    for chani = 1:size(phaDiff, 1)

        for chanj = 1:size(phaDiff, 2)

            phaDiffTr2C = squeeze(phaDiff(chani, chanj, :, :));
            meanPHDiff = angle(mean(exp(1i*(phaDiffTr2C)), 2)); %PLV

            for triali = 1:size(phaDiff, 4)

                phaDiffTr2CH = phaDiffTr2C(:, triali); 
                phaDiffTr = angle(exp(1i*(phaDiffTr2CH))); %PLV

                %disT2M{subji}(chani, chanj, :, triali) = phaDiffTr - meanPHDiff; 
                disT2M{subji}(chani, chanj, triali) = mean(phaDiffTr(2400:3200)) - mean(meanPHDiff(2400:3200)); 

            end
        end
    end
    
    
end

%% 

clear disT2
for subji = 1:10

    disT2{subji,:} = squeeze(mean(mean(disT2M{subji}))); 

end


%% VVS PFC

paths = load_paths_WM('vvs'); 
%Network_ROI_period_layers_freqs_avRep_avTimeFV_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
load ([paths.results.trial_level 'RNN_pfc_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1.mat']); 
load ([paths.electrodes_path 'pfc_elec']); 
pfc_fits = nnFit([2 3  5  9 10 11 12 14 15 16]);
pfc = pfc_fits; 



%% 
pfc_fits = pfc; 
disT2H = disT2;


tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp([2 3  5  9 10 11 12 14 15 16],1) > 1); 
%sub2exc = []; 
pfc_fits(s2e_pfc) = [];
disT2H(s2e_pfc) = []; 


ff1 = 1:6; 
tt1 = 8:13; 


nLays = 7;
clear cR

for subji = 1:length(pfc_fits)
   
    x = pfc_fits{subji};
    
    for layi = 1:7 %5:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'), 3, 'omitnan'))'; 
        y1 = disT2H{subji}; 
        cR(subji, layi, :) = corr(x1, y1, 'type', 's');
    end    
    

end


[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h)
t= squeeze(ts.tstat)








%% 



clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure(); %set(gcf, 'Position', [1000 1000 500 200])
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); hold on; colorbar; axis equal
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
%set(gca, 'FontSize', 22, 'clim', [-4 4])
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


%exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)




%% check phase differences with behavior

for subji = 1:length(disT2)

    ids = CON_ALL{subji, 9};
    
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0);
    
    idIL = cellfun(@(x) x(19), ids, 'un', 0); idIL = double(string(idIL));
    idCL = cellfun(@(x) x(20), ids, 'un', 0);idCL= double(string(idCL));
    id_Incorrect = idCL == 0 & idIL == 0;
    id_CorrIL = idIL == 1;
    id_CorrCL = idCL == 1;   
    id_CorrCNI = idIL == 0 & idCL == 1;
    disp(['Subj : ' num2str(subji) ' > CorrIL = ' num2str(sum(id_CorrIL))  ' > CorrCL = ' num2str(sum(id_CorrCL))  ' > Incorrect = ' ...
            num2str(sum(id_Incorrect)) ' > Corr category but not item = ' num2str(sum(id_CorrCNI)) ])
    nTrls(subji, 1) = sum(id_CorrIL);nTrls(subji, 2) = sum(id_CorrCL); nTrls(subji, 3) = sum(id_Incorrect);nTrls(subji, 4) = sum(id_CorrCNI);

   
    fit_CIL(subji, :) = mean(disT2{subji}(id_CorrIL), 'all');
    fit_CCL(subji, :) = mean(disT2{subji}(id_CorrCL), 'all');
    fit_Inc(subji, :) = mean(disT2{subji}(id_Incorrect), 'all');
    fit_CINC(subji, :) = mean(disT2{subji}(id_CorrCNI), 'all');
    
    
end


% % %  stats 
[h p ci ts] = ttest(fit_CIL, fit_CINC); 
t = squeeze(ts.tstat); 
disp(['p = ' num2str(p)])
disp(['t = ' num2str(ts.tstat)])





%% Across time frequency resolved 
clearvars 

tP = 2001:2800;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)

    clear PLV2U
    for freqi = 1:1:150
        t_vvs = squeeze(c_vvs.oneListTraces(:,tP,:));
        EEG = []; EEG.data = t_vvs; EEG.trials = size(t_vvs, 3); EEG.srate=1000; EEG.nbchan=size(t_vvs, 1); EEG.pnts=size(t_vvs,2);EEG.event=[];
        EEG         = pop_eegfiltnew (EEG, freqi,freqi+2);
        data_vvs        = squeeze(EEG.data); 
        dataHA_vvs      = angle(hilbert(data_vvs));
        t_pfc = squeeze(c_pfc.oneListTraces(:,tP,:));
        EEG = []; EEG.data = t_pfc; EEG.trials = size(t_pfc, 3); EEG.srate=1000; EEG.nbchan=size(t_pfc, 1); EEG.pnts=size(t_pfc,2);EEG.event=[];
        EEG         = pop_eegfiltnew (EEG, freqi,freqi+2);
        data_pfc        = squeeze(EEG.data); 
        for chani = 1:size(c_vvs.chanNames,1)
            for chanj = 1:size(c_pfc.chanNames,1)
                diffPha = angle(hilbert(squeeze(data_vvs(chani, :,:)))) - angle(hilbert(squeeze(data_pfc(chanj, :,:))));
                PLV2U(:, freqi, chani, chanj) = abs(mean(exp(1i*(diffPha))));
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


save([paths.results.PLV 'time_PLV_ALL_FR_' num2str(tP(1)) '-' num2str(tP(end))], 'PLV_ALL');




toc


%% process across time Freq Resolved

clear 
paths = load_paths_WM('vvs')
load ([paths.results.PLV 'time_PLV_ALL_FR_2001-2800'])
%load ([paths.results.PLV 'time_PLV_ALL_FR_2001-2800_1-29-30-54'])



%%
clear SI_TR MI_TR

for subji = 1:10

    allPLV = PLV_ALL{subji, 1};
    allIDs = PLV_ALL{subji, 2};
    ids0 = cellfun(@(x) strsplit(x), allIDs, 'un', 0);
    ids1 = cell2mat(cellfun(@(x) double(string(x(1))), ids0, 'un', 0));
    ids2 = cell2mat(cellfun(@(x) double(string(x(2))), ids0, 'un', 0));

    
    SI_TR{subji,1} = allPLV(ids1 == 7 & ids2 ~= 4,:,:,:); 
    MI_TR{subji,1} = allPLV(ids1 == 7 & ids2 == 4,:,:,:); 

end

%% plot 

d2pSI = cellfun(@(x) squeeze(mean(mean(mean(x, 1), 3), 4)), SI_TR, 'un', 0)
d2pMI = cellfun(@(x) squeeze(mean(mean(mean(x, 1), 3), 4)), MI_TR, 'un', 0)
c1 = cell2mat(d2pSI); 
c2 = cell2mat(d2pMI); 
%c1 = logit(c1)
%c2 = logit(c2)
md2pSI = mean(c1)
md2pMI = mean(c2)

[h p ci ts] = ttest(c1, c2)
 
md2pSI(md2pSI==0) = []; 
md2pMI(md2pMI==0) = []; 
h(isnan(h)) = []; 

hb = h; hb(hb==0) = nan; hb(hb==1) = 0.1; 
figure;
plot(md2pSI, 'r', 'LineWidth', 2); hold on
plot(md2pMI, 'b', 'LineWidth', 2)
plot(hb, 'LineWidth', 4)

legend({'SI' 'MI'})


















%% TIME resolved in a particular band

clearvars 

f2u = [3 8]

currentF = pwd;
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

nTimes = 5000; 
win_width = 500; 
mf = 100; 
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    %restrict time 
    c_vvs.oneListTraces = c_vvs.oneListTraces(:, 1001:6000,:);
    c_pfc.oneListTraces = c_pfc.oneListTraces(:, 1001:6000,:);
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            parfor triali = 1:size(c_vvs.oneListTraces,3)
                id2u = strsplit(c_vvs.oneListIds_c{triali});
                t_vvs = squeeze(c_vvs.oneListTraces(chani,:,triali));
                t_pfc = squeeze(c_pfc.oneListTraces(chanj,:,triali));
                EEG_vvs.data    = t_vvs;
                EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                EEG_vvs         = pop_eegfiltnew (EEG_vvs, f2u(1),f2u(2));
                data_vvs        = squeeze(EEG_vvs.data); 
                dataHA_vvs      = angle(hilbert(data_vvs));
                EEG_pfc.data    = t_pfc;
                EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                EEG_pfc         = pop_eegfiltnew (EEG_pfc, f2u(1),f2u(2));
                data_pfc        = squeeze(EEG_pfc.data); 
                dataHA_pfc      = angle(hilbert(data_pfc));
                diffPha = dataHA_vvs - dataHA_pfc;
                for timei = 1:bins 
                    timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                    diffPhaBIN = diffPha(timeBins);                                    
                    PLV2U(triali, chani, chanj, timei, :) = abs(mean(exp(1i*(diffPhaBIN))));
                end
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


cd (currentF)
save('PLV_ALL', 'PLV_ALL');



toc

%% Frequency resolved for cluster period

clearvars 

currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc';


tP = 2401:3200; 

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    %restrict time 
    
    clear PLV2U
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            for triali = 1:size(c_vvs.oneListTraces,3)
                triali
                id2u = strsplit(c_vvs.oneListIds_c{triali});
                t_vvs = squeeze(c_vvs.oneListTraces(chani,tP,triali));
                t_pfc = squeeze(c_pfc.oneListTraces(chanj,tP,triali));
                for freqi = 1:54
                    EEG_vvs.data    = t_vvs;
                    EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                    EEG_vvs         = pop_eegfiltnew (EEG_vvs, freqi,freqi);
                    data_vvs        = squeeze(EEG_vvs.data); 
                    dataHA_vvs      = angle(hilbert(data_vvs));
                    EEG_pfc.data    = t_pfc;
                    EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                    EEG_pfc         = pop_eegfiltnew (EEG_pfc, freqi,freqi);
                    data_pfc        = squeeze(EEG_pfc.data); 
                    dataHA_pfc      = angle(hilbert(data_pfc));
                    diffPha = dataHA_vvs - dataHA_pfc;
                    PLV2U(triali, chani, chanj, freqi, :) = abs(mean(exp(1i*(diffPha))));
                end
            end
        end
    end
    PLV_ALL{subji,1} = PLV2U; %pfc or vvs 
    PLV_ALL{subji,2} = c_vvs.oneListIds_c; %pfc or vvs 
end


cd (currentF)
save('PLV_ALL', 'PLV_ALL');


toc


%% GRANGER ONE TIME PERIOD

clearvars 

currentF = pwd;
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tP = 2401:3200;

order   =  50; % in ms
order_points   = order;


timewin_points = length(tP);


tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)
    %restrict time 
    c_vvs.oneListTraces = c_vvs.oneListTraces(:, tP,:);
    c_pfc.oneListTraces = c_pfc.oneListTraces(:, tP,:);

    clear y2x x2y
    
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            
            disp(['Sub: ' num2str(subji) ' // Chani: ' num2str(chani) ' // Chanj: ' num2str(chanj) ])
            
            clear tempdata
            tempdata(1,:,:) = c_vvs.oneListTraces(chani, :,:);
            tempdata(2,:,:) = c_pfc.oneListTraces(chanj, :,:);
            nTrials = size(tempdata,3);
            
            for triali = 1:nTrials
                tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
                tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
            end
            tempdata = reshape(tempdata,2,length(tP)*nTrials);

            % fit AR models (model estimation from bsmart toolbox)
            [Ax,Ex] = armorf(tempdata(1,:),nTrials,timewin_points,order_points);
            [Ay,Ey] = armorf(tempdata(2,:),nTrials,timewin_points,order_points);
            [Axy,E] = armorf(tempdata     ,nTrials,timewin_points,order_points);
            
            % time-domain causal estimate
            x2y(chani, chanj, :)=log(Ey/E(2,2));
            y2x(chani, chanj, :)=log(Ex/E(1,1));
                            
        
        end
    end
    

    GC{subji,1} = x2y; 
    GC{subji,2} = y2x; 
    GC{subji,3} = c_vvs.oneListIds_c; 
end


save([paths.results.PLV 'GC' '_' num2str(tP(1)) '-' num2str(tP(end)) ], 'GC');



toc

%% Process GRANGER PERIOD

%clear 
%paths = load_paths_WM('vvs')
%load ([paths.results.PLV 'GC'])
GC_V2P = GC(:, 1);
GC_P2V = GC(:, 2);


d2pP2V = cell2mat(cellfun(@(x) squeeze(mean(x, 'all')), GC_P2V, 'un', 0)')';
d2pV2P = cell2mat(cellfun(@(x) squeeze(mean(x, 'all')), GC_V2P, 'un', 0)')';
    

%% 2Bar 
data.data = [d2pP2V d2pV2P]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
%set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.32 .45] );
set(gca,'XTick',[1 2],'XTickLabel',{'', ''},'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.25 .25] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);





%% 

figure()
histogram (d2pP2V, 14)

figure()
histogram (d2pV2P, 14)



%%
GC_V2P = GC(:, 1);
GC_P2V = GC(:, 2);


d2pP2V = cell2mat(cellfun(@(x) squeeze(x), GC_P2V, 'un', 0)');
d2pV2P = cell2mat(cellfun(@(x) squeeze(x), GC_V2P, 'un', 0)');
d2pP2V = mean(d2pP2V, 2);
d2pV2P = mean(d2pV2P, 2);


%% Plot

figure;
plot(d2pP2V, 'r', 'LineWidth', 2); hold on
plot(d2pV2P, 'b', 'LineWidth', 2); hold on

set(gca, 'FontSize', 14)

%% check for the time period of network fit

c1R = mean(c1(:, 2400:3100), 2)
c2R = mean(c2(:, 2400:3100), 2)
[h p ci ts] = ttest(c1R, c2R);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);



%% GRANGER OVER TIME

clearvars 

currentF = pwd;
paths = load_paths_WM('vvs'); vvs_link = paths.traces;
paths = load_paths_WM('pfc'); pfc_link = paths.traces;

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];


order   =  30; % in ms
order_points   = order;

nTimes = 5000; 
win_width = 500; 
mf = 100; 
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
timewin_points = win_width;


tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    cd(paths.github)    

    %restrict time 
    c_vvs.oneListTraces = c_vvs.oneListTraces(:, 1001:6000,:);
    c_pfc.oneListTraces = c_pfc.oneListTraces(:, 1001:6000,:);

    clear y2x x2y
    
    for chani = 1:size(c_vvs.chanNames,1)
        for chanj = 1:size(c_pfc.chanNames,1)
            
            for timei = 1:bins 
                disp(['Sub: ' num2str(subji) ' // Chani: ' num2str(chani) ' // Chanj: ' num2str(chanj) ' // timei: ' num2str(timei)])
                timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                
                clear tempdata
                tempdata(1,:,:) = c_vvs.oneListTraces(chani, timeBins,:);
                tempdata(2,:,:) = c_pfc.oneListTraces(chanj, timeBins,:);
                nTrials = size(tempdata,3);
                
                for triali = 1:nTrials
                    tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
                    tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
                end
                tempdata = reshape(tempdata,2,length(timeBins)*nTrials);

                % fit AR models (model estimation from bsmart toolbox)
                [Ax,Ex] = armorf(tempdata(1,:),nTrials,timewin_points,order_points);
                [Ay,Ey] = armorf(tempdata(2,:),nTrials,timewin_points,order_points);
                [Axy,E] = armorf(tempdata     ,nTrials,timewin_points,order_points);
                
                % time-domain causal estimate
                x2y(chani, chanj, timei)=log(Ey/E(2,2));
                y2x(chani, chanj, timei)=log(Ex/E(1,1));
                                
            end
        end
    end
    

    GC{subji,1} = x2y; 
    GC{subji,2} = y2x; 
    GC{subji,3} = c_vvs.oneListIds_c; 
end


cd (currentF)
save([paths.results.PLV 'GC' '_time_resolved'  ], 'GC');




toc


%% 

d1p = squeeze(x2y(1, 1, :)); 

plot(d1p)

%% Process GRANGER

%clear 
%paths = load_paths_WM('vvs')
%load ([paths.results.PLV 'GC'])
GC_V2P = GC(:, 1);
GC_P2V = GC(:, 2);


d2pP2V = cell2mat(cellfun(@(x) squeeze(mean(mean(x))), GC_P2V, 'un', 0)');
d2pV2P = cell2mat(cellfun(@(x) squeeze(mean(mean(x))), GC_V2P, 'un', 0)');
    
%%
GC_V2P = GC(:, 1);
GC_P2V = GC(:, 2);


d2pP2V = cell2mat(cellfun(@(x) squeeze(x), GC_P2V, 'un', 0)');
d2pV2P = cell2mat(cellfun(@(x) squeeze(x), GC_V2P, 'un', 0)');
d2pP2V = mean(d2pP2V, 2);
d2pV2P = mean(d2pV2P, 2);


%% Plot

figure;
plot(d2pP2V, 'r', 'LineWidth', 2); hold on
plot(d2pV2P, 'b', 'LineWidth', 2); hold on

set(gca, 'FontSize', 14)

%% check for the time period of network fit

c1R = mean(c1(:, 2400:3100), 2)
c2R = mean(c2(:, 2400:3100), 2)
[h p ci ts] = ttest(c1R, c2R);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);





%%











%%