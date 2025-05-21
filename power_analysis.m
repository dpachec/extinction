%% % % % % % % % FIRST LOAD FILE - > START HERE
%%
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 

load ([paths.results.power file2load]); 

%%count channels 




clear nChans
%for subji = 1:30 %length(ALLEEG)
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if~isempty(EEG)
        nChans(subji,:) = length(EEG.chanlocs);
    end
end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])

nChans2 = nChans(any(nChans,2),:);
disp (['Elec per subject in amy: ' num2str(mean(nChans2)) ' ± ' num2str(std(nChans2))])


%% PLOT grand average for each condition
clc
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
        ids1 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ; % Cs+Cs+
        ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; % Cs+Cs- & Cs-Cs-
        
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; % Cs+Cs-
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-


        % % %   % % Acquisition only late trials
        % allACT = sum(double(string(Ev2(:, 2))) == 1); % always 72
        % t2t = allACT/2;
        % ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) > t2t;
        % ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) > t2t;

        % % %   % % Extinction only early trials
        % t2t = logical(zeros(192, 1)); 
        % t2t(73:110) = true; 
        % ids1 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) & t2t;
        % ids2 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) & t2t;

        % % %   % % Extinction only late trials
        % t2t = logical(zeros(192, 1)); 
        % t2t(111:144) = true; 
        % ids1 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) & t2t;
        % ids2 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) & t2t;

        
        % % % % Acquisition vs extinction
        %ids1 = strcmp(Ev2(:, 2), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') ;

        % % % % Acquisition vs Extinction CS+
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids2 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '1') ) ;

        % % % % Acquisition vs Extinction CS-
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '3')) ;
        %ids2 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ;

        % % % % Acquisition vs Extinction CS++
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1') ) ;
        %ids2 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '1') ) ;

        % % % % Acquisition vs Extinction CS+-
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '2') ) ;
        %ids2 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '2') ) ;

        % % % % Acquisition vs Extinction CS--
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '3') ) ;
        %ids2 = strcmp(Ev2(:, 2), '2') &  ( strcmp(Ev2(:, 6), '3') ) ;

        % % % % % % Contingencies change during EXTINCTON
        %ids1 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '3') ); % Cs++ or Cs--
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  ) ; % Cs+- 

        % % % % % % Previous plus
        %ids1 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ); % Cs++ or Cs+-
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3')  ) ; % Cs-- 

        % % % % % % CS+- vs. CS-- during extinction
        %ids1 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ); % 
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3')  ) ; % 


        nTrials(subji, :) = get_number_of_trials_POW_EXT(EEG, ids1, ids2, []); 

        



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

%freq2u = 1:12;
freq2u = 1:44;
%freq2u = 4:8;

plotC = 2; % 1: only biggest cluster; 2 all clusters; 3 no clusters

title1 = 'CS+';
title2 = 'CS-';
title3 = 'CS+ vs CS-';

%sub2excApriori = [37]'; 
%sub2excApriori = [17]'; 
sub2excApriori = []'; 

minTr = 8; 
f2u = strsplit(file2load, '_'); 


sub2exc = find(nTrials(:, 1) < minTr | nTrials(:, 2) < minTr); 
sub2exc2 = find(nTrials(:, 1) <minTr | nTrials(:, 2) <minTr ); 
sub2exc3 = [sub2excApriori; sub2exc2]; 
sub2exc = unique(sub2exc3)'; 


c1B = c1(:, freq2u, 276:475); c2B = c2(:, freq2u, 276:475); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));


[h p ci ts] = ttest(c1B, c2B, 'Alpha',0.01); 
h = squeeze(h); t = squeeze(ts.tstat);
%h(:, 1:25) = 0; 


clear allSTOBS  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTOBS(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTOBS));
max_clust_obs = allSTOBS(id); 

%

if plotC == 1
    h = zeros(size(c1B, 2),size(c1B, 3));
    h(clustinfo.PixelIdxList{id}) = 1; 
elseif plotC == 2
    %h = zeros(size(c1B, 2),size(c1B, 3));
    %h(clustinfo.PixelIdxList{id}) = 1; 
elseif plotC == 3 
    h = zeros(size(c1B, 2),size(c1B, 3));
end


hTestTimes = sum(h)'; 
hTestFreqs = sum(h, 2); 



times = -.24:.01:1.75; 
hTestTimes(:, 2) = times'; 


freqs = 1:size(c1B, 2);
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
ax1 = nexttile;
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power'; 
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125])
plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
title(title1)

xlabel('Time (s)')
ylabel('Frequency (Hz)')
ax2 = nexttile;
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);

title(title2)
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ax3 = nexttile;
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T'; c.Label.Rotation= 360; 
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
%contour(times, freqs,h, 1,':', 'Color', [0, 0, 0], 'LineWidth', 3); set(gca, 'clim', [-4 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
title(title3) 
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(ax1, brewermap([],'*Spectral'))
colormap(ax2, brewermap([],'*Spectral'))
colormap(ax3, brewermap([],'*RdYlBu'))

set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [1  length(freq2u)], 'yticklabels', {num2str(freq2u(1)), num2str(freq2u(end))}, 'xlim', [-.24 1.75]);

if freq2u(end) == 44
    set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [1  length(freq2u)], 'yticklabels', {num2str(freq2u(1)), '100'}, 'xlim', [-.24 1.75]);
end

%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'}, 'xlim', [-.25 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)
exportgraphics(gcf, ['myP.png'], 'Resolution',300)

%%nTrials2 = nTrials(any(nTrials,2),:);

%% mean in cluster 
paths = load_paths_EXT; 
load ([paths.results.additional_results 'clustInfoAMY_1-12Hz.mat'])
for subji = 1:size(c1B, 1)
    d2p1B   = squeeze(c1B(subji, :, :));
    d2p1	= mean(d2p1B(clustinfo.PixelIdxList{4}), 'all');
    d2p2B   = squeeze(c2B(subji, :, :));
    d2p2	= mean(d2p2B(clustinfo.PixelIdxList{4}), 'all');

    both2usdiff(subji, 1) = d2p1;
    both2usdiff(subji, 2) = d2p2; 
    
end



%% plot 2 bar (perform separately for Paris and Guangzhou) 

data = both2usdiff; 

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data, 1, 'omitnan');
hb = plot ([1 2], data, 'k'); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [-.25 .45] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data(:,1), data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

exportgraphics(gcf, ['myP.png'], 'Resolution',300)

%%
clc
[h p ci t] = ttest (data(:,1))


%% permutations shuffling condition labels at the group level


clearvars -except ALLEEG ids1 ids2 max_clust_obs file2load paths freq2u sub2exc allSTOBS c1B c2B nTrials c1 c2

nPerm = 1000; 
t4p = 26:200; 

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


max_clust_obs = allSTOBS(1)
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

toc

%%
max_clust_obs = allSTOBS(3)
allAb = max_clust_sum_perm((max_clust_sum_perm) > (max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% permutations shuffling trial labels in each subject


clearvars -except ALLEEG ids1 ids2 max_clust_obs file2load paths freq2u sub2exc allSTOBS c1B c2B nTrials c1 c2


nPerm = 1000; 
t4p = 301:475; 


tic
for permi = 1:nPerm
    progress_in_console(permi)
    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji};
        if ~isempty(EEG)
            c1N = squeeze(mean(EEG.power(ids1, :, freq2u ,t4p), 2, 'omitnan'));  %mean across electrodes
            c2N = squeeze(mean(EEG.power(ids2, :, freq2u ,t4p), 2, 'omitnan')); 
            junts = cat(1, c1N, c2N); 
            realCondMapping = [zeros(1,size(c1N, 1)) ones(1, size(c2N, 1))]';
            fakeCondMapping = realCondMapping(randperm(length(realCondMapping)));
            c1P(subji, :, :) = mean(junts(fakeCondMapping ==0,:,:)); 
            c2P(subji, :, :) = mean(junts(fakeCondMapping ==1,:,:)); 
        end
    end
    
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



%max_clust_obs = allSTOBS(6)
nPerm = permi; 
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

toc



%% take in cluster only
% no need to exclude subjects, already excluded in previous blocks
clearvars -except ALLEEG file2load
paths = load_paths_EXT; 
load (['_44_allTHETAPOWER'])
%load ([paths.results.additional_results '_44_allTHETAPOWER'])
load ([paths.results.additional_results '_44_clustinfo_AMY_THETA_px32'])
px = 32; 

% load _4-8_allTHETAPOWER
% load _4-8_clustinfo_AMY_THETA_px4
% px = 4; 
% 


for subji = 1:size(c1B_A, 1)

    CSPA_S = squeeze(c1B_A(subji,:,:));
    CSMA_S = squeeze(c2B_A(subji,:,:));
    CSPPE_S = squeeze(c1B_E(subji,:,:));
    CSPME_S = squeeze(c2B_E(subji,:,:));
    CSMME_S = squeeze(c3B_E(subji,:,:));
    
    CSPA(subji,:) = mean(CSPA_S(clustinfo.PixelIdxList{px}), 'all');
    CSMA(subji,:) = mean(CSMA_S(clustinfo.PixelIdxList{px}), 'all');
    CSPPE(subji,:) = mean(CSPPE_S(clustinfo.PixelIdxList{px}), 'all');
    CSPME(subji,:) = mean(CSPME_S(clustinfo.PixelIdxList{px}), 'all');
    CSMME(subji,:) = mean(CSMME_S(clustinfo.PixelIdxList{px}), 'all');

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

%%
writematrix(d2p, 'source_data_std.xlsx')


%% compare during acquisition
clc
[h p ci ts] = ttest(CSPA, CSMA); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE, CSPME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE, CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPME, CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPPE); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSPME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);

[h p ci ts] = ttest(CSMME); 
t = ts.tstat; 
disp (['t(' num2str(size(CSPA, 1)-1) ')= ' num2str(t) ',' ' p = ' num2str(p)]);


%% compare during extinction


d4ANOVA = [CSMME CSPME CSPPE];
d4ANOVA = d4ANOVA(any(d4ANOVA ,2),:);
d4ANOVA = d4ANOVA(:);
conds = [ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]';
subID = repmat(1:32, 1, 3)'; 

d4ANOVA = [d4ANOVA conds subID]

x = RMAOV1(d4ANOVA);































%%



