%% 
%% Export power

clear, close all
paths = load_paths_EXT; 

c2u = 'C';

sROI = {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' ...
        'occipital' '-cuneus' 'lingual' 'pericalcarine'}; 

%{'Amygdala'}; 
%{'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };
%{'occipital' '-cuneus' 'lingual' 'pericalcarine'}; % '-' needed for cuneus, to not be confounded with precuneus
%{'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
%{'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole' 'orbitofrontal'};

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ...
            'p_sub19', 'p_sub20', 'p_sub21'}'; %subject 8 has differnet format (see below)}';


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

        % % % remove renewal trials 
        %id2r3 = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2'); 
        %EEG.event = EEG.event(ids & id2r3);

        %epoch data and markers
        [EEG id2check] = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
        
        EEG = remove_elec_EXT(EEG, 50); %thres channels is 1/5 of 192 = 38

        

        if ~isempty(EEG.data)
            
            EEG = extract_power_EXT(EEG, 0.01); 
            %EEG = extract_theta_power_EXT(EEG); %HILBERT
            EEG = normalize_EXT(EEG);
            nChans(subji, :) = size(EEG.power, 2);
            ALLEEG{subji,:} = EEG; 

        end
    
    end


end

sROI = char(join(sROI, '_'));
mkdir(paths.results.power)
filename = [paths.results.power 'allS_' sROI '_' c2u];
nSub = sum(cell2mat(cellfun(@(x) ~isempty(x), ALLEEG, 'un', 0)));
totalChans = sum(nChans);
save(filename, 'ALLEEG', 'nSub', 'nChans', 'totalChans', '-v7.3');

cd (paths.github)


%% load traces

%file2load = ['TR_' 'HPC' '_C']; 
file2load = ['TR_' 'AMY' '_C']; 
load ([paths.results.traces file2load]); 
disp('loaded')



%% % % % % % % % FIRST LOAD FILE - > START HERE
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 

load ([paths.results.power file2load]); 

sub2exc = []; % p_07 c_27 %Always in the 50 subject space
%sub2exc = [27 37]; % p_07 c_27 %Always in the 50 subject space
ALLEEG(sub2exc,:) = {[]}; 





%% count channels 

clear nChans
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if~isempty(EEG)
        length(EEG.chanlocs)
        nChans(subji,:) = length(EEG.chanlocs);
    end
end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])



%% SELECT CHANNELS TO PLOT
clearvars -except ALLEEG paths file2load
clc 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        ids2rem1 = []; 
        %ids2rem1 = contains(chans, 'Left')
        ids2rem1 = contains(chans, 'Right')
        %ids2rem1 = contains(chans, 'Hippocampus')
        %ids2rem = logical(ids2rem1+ids2rem2)
        ids2rem = ids2rem1; 
        if ndims(EEG.power) == 3
            tmph(:, 1, :, :)  = EEG.power; 
            EEG.power = tmph; 
        end
        EEG.chanlocs(ids2rem) = []; 
        if size(EEG.chanlocs, 2) >0 & size(EEG.chanlocs, 1) >0 
            EEG.power(:, ids2rem, :, :) = []; 
            ALLEEG1{subji,:} = EEG; 
        end

    end
end

clear ALLLEEG
ALLEEG = ALLEEG1; 

%% Plot time frequency power for all channels in one subject

EEG = ALLEEG{19}; 
size(EEG.power)
nChans = size(EEG.power, 2); 
nTrials = size(EEG.power, 1); 
for chani = 1:nChans
    for triali = 1:5 %:nTrials
        
        figure()
        imagesc(squeeze(EEG.power(triali, chani, :, :))); 

    end
end


%% PLOT grand average for each condition

clearvars -except ALLEEG paths  totalChans nChans nSub nL file2load


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    if ~isempty(EEG)
        
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??


        % % %   % % Acquisition
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        ids1 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; % Cs+Cs- & Cs-Cs-
        
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; % Cs+Cs-
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-


        % % %   % % Acquisition only late trials
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) > 42;
        %ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) > 42;


        nTrials(subji, :) = get_number_of_trials_POW_EXT(EEG, ids1, ids2); 

        



        tfDCH1 = mean(EEG.power(ids1, :, : ,:), 'omitnan'); 
        tfDTF1 = squeeze(mean(tfDCH1, 2, 'omitnan'));
        tfDCH2 = mean(EEG.power(ids2, :, : ,:), 'omitnan'); 
        tfDTF2 = squeeze(mean(tfDCH2, 2, 'omitnan'));
        if exist('ids3')
            tfDCH3 = mean(EEG.power(ids3, :, : ,:), 'omitnan'); 
            tfDTF3 = squeeze(mean(tfDCH3, 2, 'omitnan'));
        end

        % % % % just to check number of trials
        %idsE1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) > 40;
        %idsL1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) <= 40;
        %xh1(subji) = sum(idsE1);         xh2(subji) = sum(idsL1); 
        %disp([num2str(sum(idsE1)) ' // ' num2str(sum(idsL1))])

        c1(subji, :, :) = tfDTF1; 
        c2(subji, :, :) = tfDTF2; 
        if exist('ids3')
            c3(subji, :, :) = tfDTF3; 
        end

                
    end

end


cd (paths.github)


%% ALL FREQUENCIES 

freq2u = 1:12;

f2u = strsplit(file2load, '_')

if strcmp(f2u{2}, 'AMY')
    sub2exc = [19 27 37]; 
elseif strcmp(f2u{2}, 'HPC')
    sub2exc = []; % no exclusions: no subjects with less than 10 trials in HPC
elseif strcmp(f2u{2}, 'OCC')
    sub2exc = []; % no exclusions: no subjects with less than 10 trials in OCC
elseif strcmp(f2u{2}, 'PFC')
    sub2exc = [19 27 37]; 
elseif strcmp(f2u{2}, 'OFC')
    sub2exc = [27]; %no subjects to exclude in OFC
elseif strcmp(f2u{2}, 'TMP')
    %sub2exc = [6 19 21 24 25 31 32 33 38 39 41 42 46 48 50]; % all subj with less than 10
    sub2exc = [19 21 25 31 46]; % all subj with less than 5 trials
end


c1B = c1(:, freq2u, 276:475); c2B = c2(:, freq2u, 276:475); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));


[h p ci ts] = ttest(c1B, c2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 

% 

h = zeros(size(c1B, 2),size(c1B, 3));
h(clustinfo.PixelIdxList{id}) = 1; 



times = -.24:.01:1.75; 
findLAT_H = sum(h);
findLAT_H(2, :) = times; 

freqs = 1:size(c1B, 2);
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
nexttile
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125])
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
title('CS+')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
nexttile
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
title('CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
title('CS+ vs CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))

set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [1  length(freq2u)], 'yticklabels', {num2str(freq2u(1)), num2str(freq2u(end))}, 'xlim', [-.24 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'}, 'xlim', [-.25 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)
exportgraphics(gcf, ['myP.png'], 'Resolution',300)


%% permutations shuffling condition labels at the group level

nPerm = 1000; 
t4p = 26:200; 

clear max_clust_sum_perm c1BP c2BP
tic
for permi = 1:nPerm
    progress_in_console(permi)
    for subji = 1:size(c1B, 1)
        if rand>.5
           c1BP(subji, :, :) = c2B(subji, :, t4p);
           c2BP(subji, :, :) = c1B(subji, :, t4p);
        else
           c1BP(subji, :, :) = c1B(subji, :, t4p); 
           c2BP(subji, :, :) = c2B(subji, :, t4p);
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1BP, c2BP); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = max(abs((allSTs)));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

toc

%% permutations shuffling trial labels in each subject


clearvars -except ALLEEG ids1 ids2 max_clust_obs file2load paths freq2u sub2exc


nPerm = 500; 
t4p = 301:475; 


tic
for subji = 1:length(ALLEEG)
    progress_in_console(subji)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        c1 = squeeze(mean(EEG.power(ids1, :, freq2u ,t4p), 2, 'omitnan'));  %mean across electrodes
        c2 = squeeze(mean(EEG.power(ids2, :, freq2u ,t4p), 2, 'omitnan')); 
        junts = cat(1, c1, c2); 
        realCondMapping = [zeros(1,size(c1, 1)) ones(1, size(c2, 1))]';
        
        for permi = 1:nPerm
            fakeCondMapping = realCondMapping(randperm(length(realCondMapping)));
            c1P = junts(fakeCondMapping ==0,:,:); 
            c2P = junts(fakeCondMapping ==1,:,:); 
        
            powPerm(subji, permi, :, :, :) = cat(1, c1P, c2P);         
        end
    end
end

powPerm(sub2exc, :, :, :, :) = []; 

toc


%%
tic
for permi = 1:nPerm
    progress_in_console(permi)
    c1P = mean(powPerm(:, permi, 1:24, :, :), 3, 'omitnan'); 
    c2P = mean(powPerm(:, permi, 25:72, :, :), 3, 'omitnan'); 

    [hPerm p ci tsPerm] = ttest(c1P, c2P); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = max(abs((allSTs)));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
end

allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm
toc

%% 
histogram(max_clust_sum_perm, 20)



%% permutations 

nPerm = 1000; 
t4p = 26:200; 


clear max_clust_sum_perm
junts = cat(1, c1B, c2B);
junts = junts(:, :, t4p); 

nSubj =  size(c1B, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';
[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);

tic
for permi = 1:nPerm
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);
    

    cond1P = junts(fakeCondMapping == 0, :,:);
    cond2P = junts(fakeCondMapping == 1, :,:);
    
    [hPerm p ci tsPerm] = ttest(cond1P, cond2P); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = max(abs((allSTs)));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


toc


%% plot cluster depic schematic 

eE1 = zeros(54,200); 
eE1(clustinfo.PixelIdxList{id}) = 1; 


figure(); set(gcf, 'Position', [100 100 500 200])
contourf(1:2000, 1:540, myresizem(eE1, 10), 40, 'linecolor', 'none'); hold on; 
contour(1:2000, 1:540,myresizem(eE1, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 6); 
plot([250 250 ],get(gca,'ylim'), 'k:','lineWidth', 8);set(gca, 'clim', [-.125 .125]); hold on; 

set(gca, 'xTick', [], 'xTicklabel', [], 'yTick', [], 'yTicklabel', []);
colormap autumn

exportgraphics(gcf, ['myP.png'], 'Resolution',300)




%% take in cluster only
clear 
load allTHETAPOWER
load clustinfo_AMY_THETA_px4


for subji = 1:size(c1B_A, 1)

    CSPA_S = squeeze(c1B_A(subji,:,:));
    CSMA_S = squeeze(c2B_A(subji,:,:));
    CSPPE_S = squeeze(c1B_E(subji,:,:));
    CSPME_S = squeeze(c2B_E(subji,:,:));
    CSMME_S = squeeze(c3B_E(subji,:,:));
    
    CSPA(subji,:) = mean(CSPA_S(clustinfo.PixelIdxList{4}), 'all');
    CSMA(subji,:) = mean(CSMA_S(clustinfo.PixelIdxList{4}), 'all');
    CSPPE(subji,:) = mean(CSPPE_S(clustinfo.PixelIdxList{4}), 'all');
    CSPME(subji,:) = mean(CSPME_S(clustinfo.PixelIdxList{4}), 'all');
    CSMME(subji,:) = mean(CSMME_S(clustinfo.PixelIdxList{4}), 'all');

end

%% 
CSMME(isnan(CSMME)) = []; 
CSPME(isnan(CSPME)) = []; 
CSPPE(isnan(CSPPE)) = []; 
CSPA(isnan(CSPA)) = []; 
CSMA(isnan(CSMA)) = []; 


d2p = [CSMME CSPME CSPPE CSPA CSMA ]; 
d2pm = mean(d2p, 'omitnan')

mean_S = mean(d2p, 'omitnan');
h = barh (mean_S);hold on;
set(h,'FaceColor', 'flat', 'lineWidth', 2);
set(gca,'YTick',[1:5],'YTickLabel',{'', ''}, 'FontSize', 18, 'linew',2, 'xlim', [-.5 .7], 'ylim', [0 6] );
plot([0 0],get(gca,'ylim'), 'k','lineWidth', 2);

fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;
h.CData = [cmap_default(5,:); cmap_default(3,:); cmap_default(2, :);  cmap_default(5,:); cmap_default(2,:)]
hb = scatter (d2p, 1:5,35, 'k'); hold on;
%h.CData = [0 0 0; 1 0 1; 1 0 1; 0 0 1; 1 0 0]

%[h p ci t] = ttest (data.data(:,1));
%disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
set(gca, 'LineWidth', 0.5, 'FontSize', 14);
box on
exportgraphics(gcf, ['myP.png'], 'Resolution',300)


%% compare during acquisition
clc
[h p ci ts] = ttest(CSPA, CSMA); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE, CSPME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE, CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPME, CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)) ')= ' num2str(t) '  ' ' p = ' num2str(p)]);


%% compare during extinction


d4ANOVA = [CSMME CSPME CSPPE];
d4ANOVA = d4ANOVA(any(d4ANOVA ,2),:);
d4ANOVA = d4ANOVA(:);
conds = [ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]';
subID = repmat(1:32, 1, 3)'; 

d4ANOVA = [d4ANOVA conds subID]

x = RMAOV1(d4ANOVA);


%% compare during acquisition

[h1 p1 ci ts] = ttest(CSPPE, CSPME); 
t1 = ts.tstat; 

[h2 p2 ci ts] = ttest(CSPPE, CSMME); 
t2 = ts.tstat; 

disp (['t = ' num2str(t1) '  ' ' p = ' num2str(p1)]);
disp (['t = ' num2str(t2) '  ' ' p = ' num2str(p2)]);



%% EXTRACT SINGLE TRIAL POWER LEVELS

clearvars -except ALLEEG paths clustinfo nL 
close all, clc
sub2exc = []; %subj 23-35-38-39-44-45

load clustinfoE_px4

for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2');
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        %length(find(ids))
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:500), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            %cTR = squeeze(mean(powH(triali, :, 3:8, 441:460), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
            nNan{subji,:} = sum(isnan(thPow));

           allPOWAMY{subji, :} = thPow; 


    end

    

end


save ('allPOWAMY', 'allPOWAMY', 'nNan'); 
 


%% Correlate responses and power in cluster

clearvars -except ALLEEG paths clustinfo nL 
close all, clc
sub2exc = []; %subj 23-35-38-39-44-45

load clustinfoE_px4


for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2');
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        %length(find(ids))
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:500), 2));
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{5}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            
            %cTR = squeeze(mean(powH(triali, :, 3:8, 441:460), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        nNan(subji,:) = sum(isnan(thPow));

        ratings2u = double(string(Ev2(ids, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        allNTR(subji,:) = length(thPow);



        allRho(subji, :) = corr(thPow, ratings2u, 'type', 'k');
        
% % % %         
% % % %         figure()
% % % %         scatter(thPow, ratings2u, 150, 'filled');
% % % %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% % % %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% % % %         allSlopes(subji, :) = C(2);
% % % %         allIntercepts(subji, :) = C(1);
% % % %         set(gca, 'ylim', [1 5], 'xlim', [-4 4], 'Fontsize', 24)
% % % %         


    end

    

end
 
idF = allRho==0 | isnan(allRho); 
allRho(idF) = []; 
%allSlopes(idF) = []; 

[h p ci t] = ttest (allRho);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);



%% plot one bar
data.data = [allRho]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.3 .5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
set(gca, 'LineWidth', 3);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)










%% ALL FREQUENCIES + ERPs

load allERPs
sub2exc = [];
c1B = c1(:, 1:54, 251:480); c2B = c2(:, 1:54, 251:480); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));


[h p ci ts] = ttest(c1B, c2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 

% 

%h = zeros(54, 230);
%h(clustinfo.PixelIdxList{id}) = 1; 

 

c2p1 = allERPs(:,1);
c2pm = cellfun(@(x) mean(x, 1), c2p1, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm1 = cat(1, c2pm{:});
cp1m = downsample(c2pm1', 10)'; 
cp1m = cp1m(:,  251:480); 
d2p1ERP = mean(cp1m); 

c2p2 = allERPs(:, 2);
c2pm = cellfun(@(x) mean(x, 1), c2p2, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm2 = cat(1, c2pm{:});
cp2m = downsample(c2pm2', 10)'; 
cp2m = cp2m(:,  251:480); 
d2p2ERP = mean(cp2m); 


times = -.5:.01:1.79; 
%times = -3:.01:3.99
freqs = 1:size(c1B, 2);
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
nexttile
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125]); hold on; 
yyaxis right
plot(times, d2p1ERP, 'k', LineWidth=2);
set(gca, ylim=[-2 2])
title('CS+')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


nexttile
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])
yyaxis right
plot(times, d2p2ERP, 'k',LineWidth=2); 
title('CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
yyaxis right
plot(times, d2p1ERP-d2p2ERP,'k', LineWidth=2); 
title('CS+ vs CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))

set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'}, 'xlim', [-.5 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power file2load '.png'], 'Resolution',150)
exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)














%% THETA BAND - LINE PLOTS

sub2exc = [];

c1B = c1(:, 3:8, 201:500); c2B = c2(:, 3:8, 201:500); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
mc1B = squeeze(mean(c1B, 2, 'omitnan'));
mc2B = squeeze(mean(c2B, 2, 'omitnan'));
mc1B(any(isnan(mc1B), 2), :) = [];
mc2B(any(isnan(mc2B), 2), :) = [];

d2pm1	= squeeze(mean(mc1B,'omitnan'));
d2pm2	= squeeze(mean(mc2B,'omitnan'));
d2pstd1	= std(mc1B);
d2pstd2	= std(mc2B);
se1 = d2pstd1/sqrt(size(mc1B, 1))
se2 = d2pstd2/sqrt(size(mc2B, 1))

[h p ci ts] = ttest(mc1B, mc2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 

h(1:150) = 0;

hb = h; hb(h==0) = nan; hb(hb==1) = -.01; 
times = (-1:.01:1.99) + .25;

colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 

%plot(times, d2p1); hold on; 
%plot(times, d2p2); hold on; 

xlabel('Time (s)')
ylabel('Theta Power')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.25 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',150)




%% permutations 2D (line plot)

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = squeeze(mc1B(:, 51:200)); 
    c2B = squeeze(mc2B(:, 51:200));
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

%%

clear p ratings2u mcsP

ratings2u = max_clust_obs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm





%% 
figure
histogram(max_clust_sum_perm); hold on; 
scatter(max_clust_obs,0, 200, 'filled','r');
set(gca, 'FontSize', 14);
xlabel('T')
exportgraphics(gcf, [paths.results.power 'myP.png'], 'Resolution',150)





%% Extract data 4 LME


clearvars -except ALLEEG paths  
close all
load clustinfoE_px5
%pi2u = unique([clustinfo_A.PixelIdxList{3} ; clustinfo_E.PixelIdxList{3} ]);
sub2exc = [];

for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        powH = EEG.power;
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:480), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{5}), 'all');
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');
            
            %cTR = squeeze(mean(powH(triali, :, 3:8, 301:400), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        ratings2u = double(string(Ev2(:, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        trialN = double(string(Ev2(:, 1)));
        exph = double(string(Ev2(:, 2)));
        currCS = double(string(Ev2(:, 8)));
        trial_type = double(string(Ev2(:, 6)));

        trialN(ids2rem) = []; 
        trial_type(ids2rem) = []; 
        currCS(ids2rem) = []; 
        exph(ids2rem) = []; 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        subID = repelem(subji, length(ratings2u))';
        
        
        d4LME = [thPow ratings2u subID trialN exph currCS trial_type];
        % % remove the testing phase for now
        %d4LME(d4LME(:, 5) == 3,:) = []; 

        allData4LME{subji,:} = d4LME;
        
    end
end


%% 

%%fit model
clc

d4LME = cat(1, allData4LME{:});
d4LME(d4LME(:, 2) == 5,:) = []; % % remove ratings of 5 !
d4LME(d4LME(:, 5) == 3,:) = []; 

tbl = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), d4LME(:,5), d4LME(:,6), d4LME(:,7),...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'Phase', 'currCS', 'trial_type'});

%tbl.Ratings = ordinal(tbl.Ratings);
%tbl2.currCS= categorical(tbl2.currCS);
%tbl.trial_type= categorical(tbl.trial_type);
%tbl2.Phase= categorical(tbl2.Phase);


%VarDecompTbl = colldiag(d4LME)
% That can be passed along for visualization
%colldiag_tableplot(VarDecompTbl);



%lme = fitlme(tbl2,'Ratings ~ theta_AMY + trial_type + Phase + trialN + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'Ratings ~ theta_AMY + currCS+ Phase+ trialN + currCS*theta_AMY + Phase*theta_AMY + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'theta_AMY ~ Ratings + currCS + Phase+ trialN + (1|subID)'); % random intercept model

lme = fitlme(tbl,'theta_AMY ~ Ratings  + trial_type + Phase*currCS + (1|subID)'); % random intercept model

lme
%lme.Coefficients


%%
writetable(tbl, 'd4LME.csv')
writetable(d4LME2, 'd4LME2.csv')


%% 
plotResiduals(lme,'fitted', 'ResidualType', 'Raw')

%plotResiduals(lme,'lagged')



%% FIT MODEL ONLY during extinction
clc

d4LME2 = cat(1, allData4LME{:});

% % % EXTINCTION
d4LME2 = d4LME2(d4LME2(:, 5) == 2, :); 
tbl2 = table(d4LME2(:,1), d4LME2(:,2), d4LME2(:,3), d4LME2(:,4),d4LME2(:,6), d4LME2(:,7), ...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'currCS', 'trial_type'});
lme = fitlme(tbl2,'theta_AMY ~ Ratings + trialN  + currCS + currCS*Ratings + (1|subID)'); % random intercept model

% % % % ACQUISITION
%d4LME2 = d4LME2(d4LME2(:, 5) == 1 & d4LME2(:, 6) == 0, :); % only cs-
%tbl2 = table(d4LME2(:,1), d4LME2(:,2), d4LME2(:,3), d4LME2(:,4), 'VariableNames',{'theta_AMY','Ratings','subID', 'trialN'});
%lme = fitlme(tbl2,'theta_AMY ~ Ratings + trialN  + (1|subID)'); % random intercept model


lme







%% 
writetable(tbl2, 'extinction_data', 'FileType','spreadsheet')


%% identify learners and no-learners

clearvars -except ALLEEG ALLEEG1 paths clustinfo

for subji = 1:size(ALLEEG1, 1)
    EEG = ALLEEG1{subji}; 
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 

        % % %   % % Acquisition
        ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        %ids1 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 

        % % %   % % both
        %ids1 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '2') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '2') ;

        mcsP(subji,:) = mean(double(string(Ev2(ids1, 7))), 'omitnan');
        mcsM(subji, :) = mean(double(string(Ev2(ids2, 7))), 'omitnan');


    end
end




%% plot 2 bar
mcsP(mcsP==0) = nan; mcsM(mcsM==0) = nan;
data.data = [mcsP mcsM ];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [0 5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

%export_fig(2, '_2.png','-transparent', '-r80');
%close all;   


%non-learners
nL = mcsP > mcsM; 







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







%% Searchlight analysis 






























%% 
d2p = squeeze(EEG.data(1, :, 1)); 
figure;plot(d2p)



%%
data2check = [EEG.data(1, :); EEG.marker_artifacts(1,:)*1000]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000);










%% plot example trial in one subject (WITH MORE THAN 1 electrode)

EEG = ALLEEG{37}; 
tr = 20; 
ch = 1; 

% figure
% d2p	= squeeze(EEG.dsPower(tr, ch, : ,:));
% myCmap = colormap(brewermap([],'YlOrRd'));
% colormap(myCmap)
% contourf(1:300, 1:54, d2p, 40, 'linecolor', 'none'); colorbar

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
d2p	= squeeze(EEG.power(tr, ch, : ,:));
contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar



%% plot two different conditions (average in all amygdala channels and across trials)


Ev = [{EEG.event.type}]';
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
Ev2 = cat(1, Ev1{:})
ids = strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1') % CS+CS+ during acquisition
d2p1	= squeeze(mean(mean(EEG.power(ids, :, : ,:), 'omitnan'), 'omitnan'));
ids = strcmp(Ev2(:, 6), '2')  & strcmp(Ev2(:, 2), '1') % CS+CS- during acquisition
d2p2	= squeeze(mean(mean(EEG.power(ids, :, : ,:), 'omitnan'), 'omitnan'));

myCmap = colormap(brewermap([],'YlOrRd')); colormap(myCmap)

tiledlayout(2, 1,'TileSpacing','compact');
nexttile
contourf(1:700, 1:54, d2p1, 40, 'linecolor', 'none'); hold on; %colorbar
nexttile
contourf(1:700, 1:54, d2p2, 40, 'linecolor', 'none'); hold on; %colorbar



%% plot all electrodes and all trials 

clear, clc

paths = load_paths_EXT; 


files2load = {['allS_' 'PFC' '_C']; ['allS_' 'AMY' '_C']; ['allS_' 'TMP' '_C']; ['allS_' 'HPC' '_C']; ['allS_' 'OCC' '_C']; ['allS_' 'OFC' '_C']}; 


for listi = 1:length(files2load)

    clearvars -except listi files2load paths

    file2load = files2load{listi}; 

    load ([paths.results.power file2load]); 

    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji}; 
        if ~isempty(EEG)
            nChans = size(EEG.chanlocs,1); 
            nTrials = size(EEG.power, 1); 
        
            for chani = 1:nChans
                for triali = 1:nTrials
                    figure
                    myCmap = colormap(brewermap([],'YlOrRd'));
                    colormap(myCmap)
                    d2p	= squeeze(EEG.power(triali, chani, : ,:));
                    contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar
                    fname = [paths.results.allSpectro '/' file2load(6:8) '/' file2load(6:end) '_' num2str(subji, '%02.f') '_' num2str(chani, '%02.f') '_' num2str(triali, '%03.f') '.png']; 
                    exportgraphics(gcf, fname, 'Resolution',50)
                    close all; 
                end
            end
        end
    end
end































%%



