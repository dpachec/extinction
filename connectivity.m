%% load traces
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C_6_4']; 
load ([paths.results.traces file2load]); 

%% Connectivity
clc
clearvars -except ALLEEG

foi = [3:8];
toi = [3001:4750]; 
toiB = [1251:3000]; 

for subji = 1:length(ALLEEG)
    disp(['Subji: ' num2str(subji)])

    EEG = ALLEEG{subji};

    clear ids2comp
    if ~isempty(EEG)
        Ev = extract_event_EXT(EEG);
        contrast = 'AMY-HPC'; 
        ids2comp = extract_ids_EXT(EEG, contrast);
        
        if ~isempty(ids2comp) 
            CON{subji, 1} = compute_connectivity_EXT(EEG.data, toi, toiB, foi, ids2comp, 'T'); % only time in this band
            CON{subji, 2} = Ev; 
        end
   end
        
    
end

save(['CON_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end))], 'CON');

%% without baseline
clc
clearvars -except CON ALLEEG
CON(any(cellfun('isempty', CON), 2), :) = [];

met2u = 4;

conH = CON; 
sub2exc = []; 
conH(sub2exc,:) = []; 

for subji = 1:length(conH)

    allC = conH{subji}(met2u,:, :); 
    allCB = conH{subji}(met2u+4,:,:); 
    Ev = conH{subji, 2}; 

    c1 = squeeze(mean(allC,3, 'omitnan')); 
    ids = find(all(c1 == 0,2));
    c1(ids, :) = []; 
    Ev (ids,:) = [];
    
    plvA(subji,:) = mean(c1(Ev(:, 2) == 1), 'all', 'omitnan') ; 
    plvE(subji,:) = mean(c1(Ev(:, 2) == 2), 'all', 'omitnan') ; 

    plvCSPA(subji,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 1), 'all', 'omitnan') ; 
    plvCSMA(subji,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 1), 'all', 'omitnan') ; 
    plvCSPE(subji,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 2), 'all', 'omitnan') ; 
    plvCSME(subji,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 2), 'all', 'omitnan') ; 
    
    plvCSPPE(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 1), 'all', 'omitnan'); 
    plvCSPME(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 2), 'all', 'omitnan'); 
    plvCSMME(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 3), 'all', 'omitnan'); 

    plvCSPPT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 1), 'all', 'omitnan'); 
    plvCSPMT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 2), 'all', 'omitnan'); 
    plvCSMMT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 3), 'all', 'omitnan'); 


end


%[h p ci ts] = ttest(plvA, plvE);
%[h p ci ts] = ttest(plvCSPA, plvCSMA);
%[h p ci ts] = ttest(plvCSPE, plvCSME);
%[h p ci ts] = ttest(plvCSPPE, plvCSMME);
%[h p ci ts] = ttest(plvCSPPT, plvCSMMT);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%% minus baseline
clc
clearvars -except CON ALLEEG
CON(any(cellfun('isempty', CON), 2), :) = [];

met2u = 1;

conH = CON; 
sub2exc = [16]; 
conH(sub2exc,:) = []; 


for subji = 1:length(conH)

   allC = conH{subji}(met2u,:, :); 
    allCB = conH{subji}(met2u+4,:,:); 
    Ev = conH{subji, 2}; 

    c1 = squeeze(mean(allC,3, 'omitnan')); 
    c1B = squeeze(mean(allCB,3, 'omitnan')); 
    ids = find(all(c1 == 0,2));
    c1(ids, :) = []; 
    c1B(ids, :) = []; 
    Ev (ids,:) = [];
    

    plvA(subji,:) = mean(c1(Ev(:, 2) == 1), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 1), 'all', 'omitnan'); 
    plvE(subji,:) = mean(c1(Ev(:, 2) == 2), 'all', 'omitnan') -  mean(c1B(Ev(:, 2) == 2), 'all', 'omitnan'); 

    plvCSPA(subji,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 1), 'all', 'omitnan') - mean(c1B(Ev(:, 8) == 1 & Ev(:, 2) == 1), 'all', 'omitnan'); 
    plvCSMA(subji,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 1), 'all', 'omitnan') - mean(c1B(Ev(:, 8) == 0 & Ev(:, 2) == 1), 'all', 'omitnan'); 
    plvCSPE(subji,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 2), 'all', 'omitnan') - mean(c1B(Ev(:, 8) == 1 & Ev(:, 2) == 2), 'all', 'omitnan'); 
    plvCSME(subji,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 2), 'all', 'omitnan') - mean(c1B(Ev(:, 8) == 0 & Ev(:, 2) == 2), 'all', 'omitnan'); 
    
    plvCSPPE(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 1), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 2 & Ev(:, 6) == 1), 'all', 'omitnan'); 
    plvCSPME(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 2), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 2 & Ev(:, 6) == 2), 'all', 'omitnan'); 
    plvCSMME(subji,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 3), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 2 & Ev(:, 6) == 3), 'all', 'omitnan'); 

    plvCSPPT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 1), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 3 & Ev(:, 6) == 1), 'all', 'omitnan'); 
    plvCSPMT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 2), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 3 & Ev(:, 6) == 2), 'all', 'omitnan'); 
    plvCSMMT(subji,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 3), 'all', 'omitnan') - mean(c1B(Ev(:, 2) == 3 & Ev(:, 6) == 3), 'all', 'omitnan'); 


end


%[h p ci ts] = ttest(plvA, plvE);
%[h p ci ts] = ttest(plvCSPA, plvCSMA);
[h p ci ts] = ttest(plvCSPE, plvCSME);
%[h p ci ts] = ttest(plvCSPPE, plvCSMME);
%[h p ci ts] = ttest(plvCSPPT, plvCSMMT);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);



%% plot one bar

ylim = [-.1 .1];
xlim = [0 3];
 
%data.data = [plvA plvE]; 
data.data = [plvCSPE plvCSME]; 
%data.data = [plvCSPPE plvCSMME]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1 2], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
%set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'   '},     'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);
exportgraphics(gcf, 'myP.png', 'Resolution',150)



%% load traces
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C_6_4']; 
load ([paths.results.traces file2load]); 

%% Connectivity FREQUENCY RESOLVED
clc
clearvars -except ALLEEG


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};

    disp(['Subji: ' num2str(subji)])
    if ~isempty(EEG)
       
        Ev = extract_event_EXT(EEG);
        ids2comp = extract_ids_EXT(EEG);
       
       if ~isempty(ids2comp) 
            CON{subji, 1} = compute_connectivity_EXT(EEG.data, [], [], [], ids2comp, 'TFR');
            CON{subji, 2} = Ev; 
        end
            
    end
        
    
end

save(['CON_TFR'], 'CON');





%% PLOT
clc
clearvars -except CON
CON(any(cellfun('isempty', CON), 2), :) = [];

met2u = 4;

for subji = 1:length(CON)

    allC = CON{subji,1}(met2u,:,:,:,:); 
    allC = allC(:,:, 3:40,:,:);
    Ev = CON{subji, 2}; 
    c1 = squeeze(mean(allC,4, 'omitnan')); 
    ids = any(any(c1,3),3); ids = ids(:, 1); ids = ~ids;  
    c1(ids, :,:) = []; 
    Ev (ids,:,:) = [];
    
    plvA(subji,:,:) = mean(c1(Ev(:, 2) == 1,:,:), 1, 'omitnan') ; 
    plvE(subji,:,:) = mean(c1(Ev(:, 2) == 2,:,:), 1, 'omitnan') ; 

    plvCSPA(subji,:,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 1,:,:), 1, 'omitnan') ; 
    plvCSMA(subji,:,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 1,:,:), 1, 'omitnan') ; 
    plvCSPE(subji,:,:) = mean(c1(Ev(:, 8) == 1 & Ev(:, 2) == 2,:,:), 1, 'omitnan') ; 
    plvCSME(subji,:,:) = mean(c1(Ev(:, 8) == 0 & Ev(:, 2) == 2,:,:), 1, 'omitnan') ; 
    
    plvCSPPE(subji,:,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 1,:,:), 1, 'omitnan') ; 
    plvCSPME(subji,:,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 2,:,:), 1, 'omitnan') ; 
    plvCSMME(subji,:,:) = mean(c1(Ev(:, 2) == 2 & Ev(:, 6) == 3,:,:), 1, 'omitnan') ; 

    plvCSPPT(subji,:,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 1,:,:), 1, 'omitnan') ; 
    plvCSPMT(subji,:,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 2,:,:), 1, 'omitnan') ; 
    plvCSMMT(subji,:,:) = mean(c1(Ev(:, 2) == 3 & Ev(:, 6) == 3,:,:), 1, 'omitnan') ; 




end




%%
c1 = plvCSPPT; 
c2 = plvCSPMT; 

sub2exc = []; %check sub17 and sub24 carefully 
% % remove baseline
c1 = c1(:, :, 1:19);
c2 = c2(:, :, 1:19);
c1B = c1 - mean(c1(:, :, 1:2), 3); 
c2B = c2 - mean(c2(:, :, 1:2), 3); 

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

h = zeros(38,19);
%h(clustinfo.PixelIdxList{id}) = 1; 


times = 1:190;
freqs = 1:size(c1B, 2)*10;
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
nexttile
imagesc(times, freqs, flipud(myresizem(d2p1, 10))); hold on; c = colorbar; %c.Label.String = 'Z-Power';
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.05 .05])

%contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
%plot([1.5 1.5 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125])

nexttile
imagesc(times, freqs, flipud(myresizem(d2p2, 10))); hold on; c = colorbar; %c.Label.String = 'Z-Power';
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.05 .05])

%contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
%plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])

nexttile

imagesc(times, freqs, flipud(myresizem(t, 10))); hold on; c = colorbar; %c.Label.String = 'Z-Power';
contour(times, freqs,flipud(myresizem(h,10)), 1, 'Color', [0, 0, 0], 'LineWidth', 2); 
plot([25 25],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-4 4])

%contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
%contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
%plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);

colormap(brewermap([],'*Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 380], 'yticklabels', {'40', '3'}, 'xtick', [25 75 125 175], 'xticklabels', {'0', '0.5', '1', '1.5'});

exportgraphics(gcf, [ 'myP.png'], 'Resolution',150)



%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1C = c1B(:,:,3:end); 
    c2C = c2B(:,:,3:end); 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1C(subji, :, :);
           c1C(subji, :, :) = c2C(subji, :, :);
           c2C(subji, :, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1C, c2C); 
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
ratings2u = max_clust_obs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm




%% GC analysis > load traces
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C_6_4']; 
load ([paths.results.traces file2load]); 


%% GRANGER Temporally resolved

clearvars -except ALLEEG

nTimepoints = 2500; %%all2all file is stored within this cell array
win_width   = 500; 
mf          = 50; 
bins        =  floor ( (nTimepoints/mf)- win_width/mf+1 );

order   =  50; % in ms
order_points   = order;


tic
for subji = 1:length(ALLEEG)
    disp(['Subji: ' num2str(subji)])
    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = extract_event_EXT(EEG);
        ids2comp = extract_ids_EXT(EEG);

         if ~isempty(ids2comp)

            %select specific condition
            CoI = Ev(:, 2) == 2 & Ev(:, 8) == 1; 
            EEG.data = EEG.data(:,2501:5000, CoI); 
            Ev = Ev(CoI);

            % % % remove nan trials 
            count = 1; 
            for triali = 1:size(EEG.data, 3)
                if find(isnan(EEG.data(:, :, triali)))
                    id2rem(count, :) = triali; 
                    count = count+1; 
                end
            end
            
            EEG.data(:, :, id2rem) = []; 
            Ev(id2rem,:) = []; 
            EEG.trials = size(EEG.data, 3);


            tf_granger=zeros(2,length(ids2comp(:, 1)),bins);
            tf_grangerCombi1 = zeros(length(ids2comp(:, 1)), bins);
            tf_grangerCombi2 = zeros(length(ids2comp(:, 1)), bins);
           
            parfor combi = 1:length(ids2comp(:, 1))

                chan1 = ids2comp(combi, 1);
                chan2 = ids2comp(combi, 2);
                %eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));
                eegdata = EEG.data([chan1 chan2],:,:);

                disp(['Comb: ' num2str(combi) ' of ' num2str(length(ids2comp(:, 1)))])
                
                y2x = zeros (1, bins)
                x2y = zeros (1, bins)
                for timei=1:bins
                    timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1; 
                    tempdata = squeeze(eegdata(:,timeBinsi,:));
                    for triali=1:size(tempdata, 3)
                        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
                        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
                    end
                    % reshape tempdata for armorf
                    tempdata = reshape(tempdata,2,win_width*EEG.trials);
                    % fit AR models (model estimation from bsmart toolbox)
                    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,win_width,order_points);
                    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,win_width,order_points);
                    [Axy,E] = armorf(tempdata     ,EEG.trials,win_width,order_points);
                    
                    % time-domain causal estimate
                    y2x(timei)=log(Ex/E(1,1));
                    x2y(timei)=log(Ey/E(2,2));
                    
                end
                tf_grangerCombi1(combi, :) = y2x; 
                tf_grangerCombi2(combi, :) = x2y; 
        end
        tf_granger(1,:,:)  = tf_grangerCombi1;
        tf_granger(2,:,:)  = tf_grangerCombi2;

        GC{subji,1} = tf_granger; 
        GC{subji,2} = Ev; 
    end
    end
end


save(['GCT' ], 'GC');



toc

%%
clc
clearvars -except GC
GC(any(cellfun('isempty', GC), 2), :) = [];


sub2exc = []; 


for subji = 1:length(GC)
    m1(subji,:)= squeeze(mean(GC{subji}(1,:,:), 2)); 
    m2(subji,:)= squeeze(mean(GC{subji}(2,:,:), 2)); 
end 
m1(sub2exc,:) = []; m2(sub2exc,:) = []; 

d2pm1	= squeeze(mean(m1,'omitnan'));
d2pm2	= squeeze(mean(m2,'omitnan'));
d2pstd1	= std(m1);
d2pstd2	= std(m2);
se1 = d2pstd1/sqrt(size(m1, 1));
se2 = d2pstd2/sqrt(size(m2, 1));

[h p ci ts] = ttest(m1, m2); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 



hb = h; hb(h==0) = nan; hb(hb==1) = -.01; 
times = 1:41;

colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot (times, hb, 'Linewidth', 7)









%% PSI analysis > load traces
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C_6_4']; 
load ([paths.results.traces file2load]); 

%% PSI FIELDTRIP
clearvars -except ALLEEG

toi = [3001:4750];
%toi = [3501:4750];


tic
for subji = 1:length(ALLEEG)
    disp(['Subji: ' num2str(subji)])
    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = extract_event_EXT(EEG);
        contrast = 'AMY-HPC';
        ids2comp = extract_ids_EXT(EEG, contrast);
       
        
         if ~isempty(ids2comp)
            nCombi = length(ids2comp(:, 1)); 
            %select specific condition
            CoI = Ev(:, 2) == 2 & Ev(:, 8) == 0; 
            EEG.data = EEG.data(:,toi, CoI); 
            Ev = Ev(CoI);

            % % % remove nan trials 
            count = 1; 
            for triali = 1:size(EEG.data, 3)
                if find(isnan(EEG.data(:, :, triali)))
                    id2rem(count, :) = triali; 
                    count = count+1; 
                end
            end            
            EEG.data(:, :, id2rem) = []; 
            Ev(id2rem,:) = []; 
            
            EEG = add_EEGLAB_fields(EEG);
            dataFT = eeglab2fieldtrip(EEG, 'preprocessing', 'raw');
           
            cfg = [];
            cfg.method = 'mtmfft';
            cfg.taper = 'hanning';
            cfg.foi = 1:1:40; % 1-40 Hz
            cfg.tapsmofrq = 1; % 1-Hz half-band
            cfg.pad = 5; % pad to 4 s
            cfg.output = 'fourier';
            cfg.keeptrials = 'yes';
            dataPSI = ft_freqanalysis(cfg, dataFT); % spectral decomposition
            
            cfgp = [];
            cfgp.method = 'psi';
            cfgp.bandwidth = 5;
            PSIC = ft_connectivityanalysis(cfgp, dataPSI);


            end
        
        PSI{subji,1} = PSIC; 
        PSI{subji,2} = Ev; 
    end
    end
end


save(['PSIT_' contrast ], 'PSI');



toc



%% PSI Temporally resolved
clearvars -except ALLEEG

toi = [3001:4750];
%toi = [3501:4750];


tic
for subji = 1:length(ALLEEG)
    disp(['Subji: ' num2str(subji)])
    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = extract_event_EXT(EEG);
        contrast = 'AMY-HPC';
        ids2comp = extract_ids_EXT(EEG, contrast);
       
        
         if ~isempty(ids2comp)
            nCombi = length(ids2comp(:, 1)); 
            %select specific condition
            CoI = Ev(:, 2) == 2 & Ev(:, 8) == 0; 
            EEG.data = EEG.data(:,toi, CoI); 
            Ev = Ev(CoI);

            % % % remove nan trials 
            count = 1; 
            for triali = 1:size(EEG.data, 3)
                if find(isnan(EEG.data(:, :, triali)))
                    id2rem(count, :) = triali; 
                    count = count+1; 
                end
            end            
            EEG.data(:, :, id2rem) = []; 
            Ev(id2rem,:) = []; 
            EEG.trials = size(EEG.data, 3);
            
            PSIC = zeros(nCombi,40);
            parfor combi = 1:nCombi

                chan1 = ids2comp(combi, 1);
                chan2 = ids2comp(combi, 2);
                eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));
                %eegdata = EEG.data([chan1 chan2],:,:);

                %disp(['Comb: ' num2str(combi) ' of ' num2str(nCombi)])
                PSIF = zeros(1,38);
                for freqi = 3:40
                    M =  data2psiX(eegdata, 1000, [freqi-2 freqi+2], 0);
                    PSIF(freqi) = M(1,2); %- M(1,2);
                end
                PSIC(combi,:) = PSIF;
            end
        
        PSI{subji,1} = PSIC; 
        PSI{subji,2} = Ev; 
    end
    end
end


save(['PSIT_' contrast ], 'PSI');



toc





%% PLOT
clc
clearvars -except PSI ALLEEG
PSI(any(cellfun('isempty', PSI), 2), :) = [];

clear c1 c2
for subji = 1:length(PSI)

    allC = PSI{subji,1}; 
    Ev = PSI{subji, 2}; 
    c1(subji,:) = squeeze(mean(allC,1, 'omitnan')); 
    
end

 
d2p1 = squeeze(mean(c1));
std1	= std(c1);
se1 = std1/sqrt(size(c1, 1));

[h p ci ts] = ttest(d2p1); 
h = squeeze(h); t = squeeze(ts.tstat);


shadedErrorBar(1:40,  d2p1, se1, 'k', 1); hold on; 

%% PLOT
clc
PSI_CSP(any(cellfun('isempty', PSI_CSP), 2), :) = [];
PSI_CSM(any(cellfun('isempty', PSI_CSM), 2), :) = [];

clear c1 c2
for subji = 1:length(PSI_CSP)

    allC = PSI_CSP{subji,1}; 
    c1(subji,:) = squeeze(mean(allC,1, 'omitnan')); 
    allC = PSI_CSM{subji,1}; 
    c2(subji,:) = squeeze(mean(allC,1, 'omitnan')); 
    
end

 
d2p1 = squeeze(mean(c1));
std1	= std(c1);
se1 = std1/sqrt(size(c1, 1));
d2p2 = squeeze(mean(c2));
std2	= std(c2);
se2 = std2/sqrt(size(c2, 1));


[h p ci ts] = ttest(c1, c2); 
h = squeeze(h); t = squeeze(ts.tstat);


shadedErrorBar(1:40,  d2p1, se1, 'k', 1); hold on; 
shadedErrorBar(1:40,  d2p2, se2, 'r', 1); hold on; 




%%
clear all;
 	
% simulate some data
fsample = 1000;
nsample = fsample*30;
 
dat = randn(1,nsample+100);
 
data.trial{1}(1,:) = dat(1,1:nsample) + 0.1.*randn(1,nsample);
data.trial{1}(2,:) = dat(1,100+(1:nsample)) + 0.1.*randn(1,nsample);
data.time{1} = (1:nsample)./fsample;
data.label   = {'a';'b'}; % channel 2 is leading channel 1
 
figure;plot(data.time{1},data.trial{1}); xlim([0 1]);
 
% cut into 2 second snippets
cfg = [];
cfg.length = 2;
data = ft_redefinetrial(cfg, data);
 
% spectral decomposition
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.tapsmofrq = 2;
cfg.foilim = [0 100];
freq = ft_freqanalysis(cfg, data);
 
% connectivity estimation
cfg = [];
cfg.method = 'psi';
cfg.bandwidth = 5;
psi = ft_connectivityanalysis(cfg, freq);
 
% visualization
cfg = [];
cfg.parameter = 'psispctrm';
ft_connectivityplot(cfg, psi);

%% GRANGER ONE TIME PERIOD

clearvars -except ALLEEG

toi = [3001:4750]; 

order   =  36; % in ms
order_points   = order;


timewin_points = length(toi);

min_freq = 3; % in Hz, using minimum of 10 Hz because of 200-ms window
max_freq = 40; 

order_points = 15;

frequencies = logspace(log10(min_freq),log10(max_freq),15);

% initialize
tf_granger=zeros(2,length(frequencies));


tic
for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = [{EEG.event.type}]'; 
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        clen = cellfun(@length, Ev1); 
        EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' ');
        Ev3 = double(string(Ev2));
        
        

        % % % remove nan trials 
        count = 1; 
        for triali = 1:size(EEG.data, 3)
            if find(isnan(EEG.data(:, :, triali)))
                id2rem(count, :) = triali; 
                count = count+1; 
            end
        end
        
        EEG.data(:, :, id2rem) = []; 
        Ev3(id2rem,:) = []; 

        chansLab = {EEG.chanlocs.fsLabel}';
        chAMYid = contains(chansLab, 'Amygdala'); 
        chHPCid = contains(chansLab, 'Hippocampus'); 
        chAMY = EEG.data(chAMYid, toi, :);
        chHPC = EEG.data(chHPCid, toi, :);

        clear y2x x2y
        for chani = 1:size(chAMY,1)
            for chanj = chani:size(chHPC,1)
                
                disp(['Sub: ' num2str(subji) ' // Chani: ' num2str(chani) ' // Chanj: ' num2str(chanj) ])
                
                clear tempdata
                tempdata(1,:,:) = chAMY(chani, :,:);
                tempdata(2,:,:) = chHPC(chanj, :,:);
                nTrials = size(tempdata,3);
                
                for triali = 1:nTrials
                    tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
                    tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
                end
                tempdata = reshape(tempdata,2,length(toi)*nTrials);
    
                % fit AR models (model estimation from bsmart toolbox)
                [Ax,Ex] = armorf(tempdata(1,:),nTrials,timewin_points,order_points);
                [Ay,Ey] = armorf(tempdata(2,:),nTrials,timewin_points,order_points);
                [Axy,E] = armorf(tempdata     ,nTrials,timewin_points,order_points);
                
                % code below is adapted from bsmart toolbox function pwcausal.m
                % corrected covariance
                eyx = E(2,2) - E(1,2)^2/E(1,1); 
                exy = E(1,1) - E(2,1)^2/E(2,2);
                N = size(E,1);
                
                for fi=1:length(frequencies)
                    
                    % transfer matrix (note the similarity to Fourier transform)
                    H = eye(N);
                    for m = 1:order_points
                        H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
                    end
                    
                    Hi = inv(H);
                    S  = H\E*Hi'/EEG.srate;
                    
                    % granger prediction per frequency
                    tf_granger(1,fi,:) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) );
                    tf_granger(2,fi,:) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) );
                end
            
            end
        end
    
        GC{subji,1} = tf_granger; 
        GC{subji,2} = Ev3; 
    end
end


save(['GC' '_' num2str(toi(1)) '-' num2str(toi(end)) ], 'GC');



toc



%%
clc
clearvars -except GC
GC(any(cellfun('isempty', GC), 2), :) = [];


for subji = 1:length(GC)

    GCh = GC{subji, 1}; 
    Ev = GC{subji,2}; 
    plvCSP(subji,:) = mean(GCh(Ev(:, 8) == 1,:,:), 'all', 'omitnan'); 
    plvCSM(subji,:) = mean(GCh(Ev(:, 8) == 0,:,:), 'all', 'omitnan'); 
    plvA(subji,:) = mean(GCh(Ev(:, 2) == 1,:,:), 'all', 'omitnan'); 
    plvE(subji,:) = mean(GCh(Ev(:, 2) == 2,:,:), 'all', 'omitnan'); 

    plvCSPA(subji,:) = mean(GCh(Ev(:, 8) == 1 & Ev(:, 2) == 1,:,:), 'all', 'omitnan'); 
    plvCSMA(subji,:) = mean(GCh(Ev(:, 8) == 0 & Ev(:, 2) == 1,:,:), 'all', 'omitnan'); 
    plvCSPE(subji,:) = mean(GCh(Ev(:, 8) == 1 & Ev(:, 2) == 2,:,:), 'all', 'omitnan'); 
    plvCSME(subji,:) = mean(GCh(Ev(:, 8) == 0 & Ev(:, 2) == 2,:,:), 'all', 'omitnan'); 
    
    plvCSPPE(subji,:) = mean(GCh(Ev(:, 2) == 2 & Ev(:, 6) == 1,:,:), 'all', 'omitnan'); 
    plvCSPME(subji,:) = mean(GCh(Ev(:, 2) == 2 & Ev(:, 6) == 2,:,:), 'all', 'omitnan'); 
    plvCSMME(subji,:) = mean(GCh(Ev(:, 2) == 2 & Ev(:, 6) == 3,:,:), 'all', 'omitnan'); 

    plvCSPPT(subji,:) = mean(GCh(Ev(:, 2) == 3 & Ev(:, 6) == 1,:,:), 'all', 'omitnan'); 
    plvCSPMT(subji,:) = mean(GCh(Ev(:, 2) == 3 & Ev(:, 6) == 2,:,:), 'all', 'omitnan'); 
    plvCSMMT(subji,:) = mean(GCh(Ev(:, 2) == 3 & Ev(:, 6) == 3,:,:), 'all', 'omitnan'); 


end

[h p ci ts] = ttest(plvCSP, plvCSM);
%[h p ci ts] = ttest(plvCSPE, plvCSME);
disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);




%% GC analysis > load traces
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C_6_4']; 
load ([paths.results.traces file2load]); 



%% ESTIMATE MODEL ORDER 

% load sample EEG data
EEG = ALLEEG{3};

% define channels for granger prediction
chan1name = 'POL A3';
chan2name = 'POL H1';

% Granger prediction parameters
timewin = 200; % in ms
order   =  27; % in ms

% convert parameters to indices
timewin_points = round(timewin/(1000/EEG.srate));
order_points   = round(order/(1000/EEG.srate));

% find the index of those channels
chan1 = find(strcmpi(chan1name,{EEG.chanlocs.labels}));
chan2 = find(strcmpi(chan2name,{EEG.chanlocs.labels}));

Ev = [{EEG.event.type}]'; 
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
clen = cellfun(@length, Ev1); 
EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
Ev2 = cat(1, Ev1{:});
Ev2(:, 10) = erase(Ev2(:, 10), ' ');
Ev3 = double(string(Ev2));

% % % remove nan trials 
count = 1; 
for triali = 1:size(EEG.data, 3)
    if find(isnan(EEG.data(:, :, triali)))
        id2rem(count, :) = triali; 
        count = count+1; 
    end
end

EEG.data(:, :, id2rem) = []; 
Ev3(id2rem,:) = []; 
EEG.trials = size(EEG.data, 3);

eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));
times2saveidx = 3001:20:4000;

[x2y,y2x] = deal(zeros(1,length(times2saveidx))); % the function deal assigns inputs to all outputs
bic = zeros(length(times2saveidx),70); % Bayes info criteria (hard-coded to order=15)

for timei=1:length(times2saveidx)
    timei
    times2sav(timei,:) = times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2);
    tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
    for triali=1:size(tempdata, 3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
    end
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
    parfor bici=1:size(bic,2)
        [Axy,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
        bic(timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
    end
end

% draw lines
figure
plot(times2saveidx,x2y)
hold on
plot(times2saveidx,y2x,'r')
legend({[ 'GP: ' chan1name ' -> ' chan2name ];[ 'GP: ' chan2name ' -> ' chan1name ]})
title([ 'Window length: ' num2str(timewin) ' ms, order: ' num2str(order) ' ms' ])
xlabel('Time (ms)')
ylabel('Granger prediction estimate')


%% Plot model order (Figure 28.4)

% see "bici" for-loop above for code to compute BIC

figure

subplot(121)
plot((1:size(bic,2))*(1000/EEG.srate),mean(bic,1),'--.')
xlabel('Order (converted to ms)')
ylabel('Mean BIC over all time points')

[bestbicVal,bestbicIdx]=min(mean(bic,1));
hold on
plot(bestbicIdx*(1000/EEG.srate),bestbicVal,'mo','markersize',15)

title([ 'Optimal order is ' num2str(bestbicIdx) ' (' num2str(bestbicIdx*(1000/EEG.srate)) ' ms)' ])

subplot(122)
[junk,bic_per_timepoint] = min(bic,[],2);
plot(times2saveidx,bic_per_timepoint*(1000/EEG.srate),'--.')
xlabel('Time (ms)')
ylabel('Optimal order (converted to ms)')
title('Optimal order (in ms) at each time point')





%% Figure 28.5 SPECTRAL GRANGER EXAMPLE

clearvars -except ALLEEG
EEG = ALLEEG{6};

% define channels for granger prediction
chan1name = 'POL A3';
chan2name = 'POL H1';

% Granger prediction parameters
timewin = 200; % in ms
order   =  30; % in ms
min_freq = 10; % in Hz, using minimum of 10 Hz because of 200-ms window
max_freq = 40;

frequencies = logspace(log10(min_freq),log10(max_freq),15);

% convert parameters to indices
timewin_points = round(timewin/(1000/EEG.srate));
order_points   = round(order/(1000/EEG.srate));

% find the index of those channels
chan1 = find(strcmpi(chan1name,{EEG.chanlocs.labels}));
chan2 = find(strcmpi(chan2name,{EEG.chanlocs.labels}));

Ev = extract_event_EXT(EEG);
order_points = order;

% % % remove nan trials 
count = 1; 
for triali = 1:size(EEG.data, 3)
    if find(isnan(EEG.data(:, :, triali)))
        id2rem(count, :) = triali; 
        count = count+1; 
    end
end

EEG.data(:, :, id2rem) = []; 
Ev(id2rem,:) = []; 
EEG.trials = size(EEG.data, 3);

eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));
times2saveidx = 3001:20:4000;

% initialize
tf_granger=zeros(2,length(frequencies),length(times2saveidx));



for timei=1:length(times2saveidx)
    timei
    % data from all trials in this time window
    % (note that the ERP-subtracted data are used)
    tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));

    % detrend and zscore all data
    for triali=1:size(tempdata,3)
        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));

        % At this point with real data, you might want to check for stationarity
        % and possibly discard or mark data epochs that are non-stationary.
    end

    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);

    % fit AR models
    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
    [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);

    % code below is adapted from bsmart toolbox function pwcausal.m
    % corrected covariance
    eyx = E(2,2) - E(1,2)^2/E(1,1);
    exy = E(1,1) - E(2,1)^2/E(2,2);
    N = size(E,1);

    for fi=1:length(frequencies)

        % transfer matrix (note the similarity to Fourier transform)
        H = eye(N);
        for m = 1:order_points
            H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
        end

        Hi = inv(H);
        S  = H\E*Hi'/EEG.srate;

        % granger prediction per frequency
        tf_granger(1,fi,timei) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) );
        tf_granger(2,fi,timei) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) );
    end
end % end time loop

%%
% - plot - %
figure, set(gcf,'Name',[ 'Granger prediction between electrodes' chan1name ' and ' chan2name '.' ]);
subplot(211)
contourf(times2saveidx,frequencies,squeeze(tf_granger(1,:,:)),40,'linecolor','none')
%set(gca,'clim',[-.01 .01])
colorbar
title([ chan1name '->' chan2name '; ''raw'' units' ])

subplot(212)
contourf(times2saveidx,frequencies,squeeze(tf_granger(2,:,:)),40,'linecolor','none')
%set(gca,'clim',[0 .5])
colorbar
title([ chan2name '->' chan1name '; ''raw'' units' ])


%% SPECTRAL GRANGER ALL SUBJECTS 

clearvars -except ALLEEG

timewin = 200; % in ms
order   =  30; % in ms
min_freq = 10; % in Hz, using minimum of 10 Hz because of 200-ms window
max_freq = 40;

frequencies = logspace(log10(min_freq),log10(max_freq),15);

timewin_points = timewin;
order_points = order;
times2saveidx = 3001:20:4000;

for subji = 1:length(ALLEEG)
    disp(['Subji: ' num2str(subji)])

    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = extract_event_EXT(EEG);
        ids2comp = extract_ids_EXT(EEG);

        if ~isempty(ids2comp)
            % % % remove nan trials 
            count = 1; 
            for triali = 1:size(EEG.data, 3)
                if find(isnan(EEG.data(:, :, triali)))
                    id2rem(count, :) = triali; 
                    count = count+1; 
                end
            end
            
            EEG.data(:, :, id2rem) = []; 
            Ev(id2rem,:) = []; 
            EEG.trials = size(EEG.data, 3);
        
            %tf_granger=zeros(2,length(frequencies),length(times2saveidx));
            tf_grangerCombi1 = zeros(length(ids2comp(:, 1)), length(frequencies),length(times2saveidx));
            tf_grangerCombi2 = zeros(length(ids2comp(:, 1)), length(frequencies),length(times2saveidx));
           parfor combi = 1:length(ids2comp(:, 1))
                
                disp(['Comb: ' num2str(combi) '/' num2str(length(ids2comp(:, 1))) ])
                chan1 = ids2comp(combi, 1);
                chan2 = ids2comp(combi, 2);
                
                %eegdata = bsxfun(@minus,EEG.data([chan1 chan2],:,:),mean(EEG.data([chan1 chan2],:,:),3));
                eegdata = EEG.data([chan1 chan2],:,:);
                tf_granger1 = zeros(length(frequencies), length(times2saveidx));
                tf_granger2  = zeros(length(frequencies), length(times2saveidx));
                for timei=1:length(times2saveidx)
                    % data from all trials in this time window % (note that the ERP-subtracted data are used)
                    tempdata = squeeze(eegdata(:,times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2),:));
                
                    % detrend and zscore all data
                    for triali=1:size(tempdata,3)
                        tempdata(1,:,triali) = zscore(detrend(squeeze(tempdata(1,:,triali))));
                        tempdata(2,:,triali) = zscore(detrend(squeeze(tempdata(2,:,triali))));
                        % At this point with real data, you might want to check for stationarity
                        % and possibly discard or mark data epochs that are non-stationary.
                    end
                
                    % reshape tempdata for armorf
                    tempdata = reshape(tempdata,2,timewin_points*EEG.trials);
                
                    % fit AR models
                    [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
                    [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
                    [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
                
                    % code below is adapted from bsmart toolbox function pwcausal.m
                    eyx = E(2,2) - E(1,2)^2/E(1,1);
                    exy = E(1,1) - E(2,1)^2/E(2,2);
                    N = size(E,1);
                    tf_granger1T = zeros(1,length(frequencies));
                    tf_granger2T = zeros(1,length(frequencies));
                    for fi=1:length(frequencies)
                        % transfer matrix (note the similarity to Fourier transform)
                        H = eye(N);
                        for m = 1:order_points
                            H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
                        end
                
                        Hi = inv(H);
                        S  = H\E*Hi'/EEG.srate;
                
                        % granger prediction per frequency
                        tf_granger1T(fi) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) );
                        tf_granger2T(fi) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) );
                    end
                    tf_granger1(:,timei) = tf_granger1T; 
                    tf_granger2(:,timei) = tf_granger2T; 
                end % end time loop

                tf_grangerCombi1(combi,:,:) = tf_granger1;
                tf_grangerCombi2(combi,:,:) = tf_granger2;
            end

            tf_granger(:, 1, :, :) = tf_grangerCombi1;
            tf_granger(:, 2, :, :) = tf_grangerCombi2;
            GC{subji, 1} = tf_granger; 
            GC{subji, 2} = Ev; 
        end
    end
end

save('GC', 'GC')






%% PLOT
clc
clearvars -except GC
GC(any(cellfun('isempty', GC), 2), :) = [];

clear c1 c2
for subji = 1:length(GC)

    allC = GC{subji,1}; 
    Ev = GC{subji, 2}; 
    c1(subji, :,:) = squeeze(mean(allC(:,1,:,:),1, 'omitnan')); 
    c2(subji, :,:) = squeeze(mean(allC(:,2,:,:),1, 'omitnan')); 
    
end

%% 
d2p1 = squeeze(mean(c1));
d2p2 = squeeze(mean(c2));

[h p ci ts] = ttest(c1, c2);
h = squeeze(h); t = squeeze(ts.tstat);

tiledlayout(3,1)
nexttile
contourf(1:50, 1:15, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
nexttile
contourf(1:50, 1:15, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
nexttile
contourf(1:50, 1:15, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
contour(1:50, 1:15,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-3 3])


%%
c1 = plvCSPE; 
c2 = plvCSME; 

sub2exc = []; %check sub17 and sub24 carefully 
% % remove baseline
c1 = c1(:, :, 1:19);
c2 = c2(:, :, 1:19);
c1B = c1 - mean(c1(:, :, 1:2), 3); 
c2B = c2 - mean(c2(:, :, 1:2), 3); 

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

h = zeros(50,19);
h(clustinfo.PixelIdxList{id}) = 1; 

times = 1:19; 
freqs = 1:size(c1B, 2);
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
nexttile
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
%plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125])
nexttile
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
%plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
%plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))

set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});

%exportgraphics(gcf, [ 'myP.png'], 'Resolution',150)



%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1C = c1B(:,:,3:end); 
    c2C = c2B(:,:,3:end); 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1C(subji, :, :);
           c1C(subji, :, :) = c2C(subji, :, :);
           c2C(subji, :, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1C, c2C); 
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
ratings2u = max_clust_obs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm
































%%



















%%