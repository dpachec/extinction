%% 
%% 
clear, close all
paths = load_paths_EXT; 

c2u = 'C';

%sROI = {'orbitofrontal'}; 

%sROI = {'lateraloccipital'}; 

%sROI = {'superiorfrontal' 'rostralmiddlefrontal' 'anteriorcingulate' 'posteriorcingulate' 'precentral' 'caudalmiddlefrontal'}; % case sensitive 

 sROI = { 'inferiortemporal' 'middletemporal' 'superiortemporal' 'bankssts' 'fusiform' 'temporalpole' ...
              'lateraloccipital' 'lingual' 'parahippocampal' 'cuneus' 'pericalcarine' };

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

    %chansLab = {EEG.chanlocs.fsLabelsR}';
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

        %epoch data and markers
        [EEG id2check] = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
        
        EEG = remove_elec_EXT(EEG, 50); %thres channels is 1/5 of 192 = 38
        

        if ~isempty(EEG.data)
            
            EEG = extract_power_EXT(EEG, 0.01); 
            EEG = normalize_EXT(EEG);
            EEG = rmfield(EEG, 'data');
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




%% plot example trial in one subject (ONLY 1 electrode)

EEG = ALLEEG{3}; 
tr =20; 

figure
d2p	= squeeze(EEG.power(tr, 1 ,:,:));
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar



%% plot Mean across trials in one subject (ONLY 1 electrode)

EEG = ALLEEG{3};

figure
d2p	= squeeze(mean(EEG.power(:, 1 ,:,:), 'omitnan'));
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar



%% delete trigger labels in events

x = [{EEG2.event.type}]; 
ids2rem = strcmp(x, 'trigger');
EEG2.event(ids2rem) = []; 

%% check that markers are ok
%data2check = [EEG.data(1, :); EEG.markers_artifacts(1,:)*1000]; 
data2check = [EEG.data(1, :)]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000, 'events', EEG.event);



%% % % % % % % % FIRST LOAD FILE - > START HERE
clear

paths = load_paths_EXT; 
file2load = ['allS_' 'Amygdala' '_C']; 
%file2load = ['allS_' 'Hippocampus' '_C']; 
%file2load = ['allS_' 'orbitofrontal' '_C']; 
load ([paths.results.power file2load]); 






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



%% SELECT CHANNELS TO PLOT
clearvars -except ALLEEG paths file2load
clc 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        ids2rem1 = []; 
        %ids2rem1 = contains(chans, 'Right')
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

%% PLOT grand average for each condition

clearvars -except ALLEEG ALLEEG1 paths  totalChans nChans nSub nL



for subji = 1:length(ALLEEG1)
    
    EEG = ALLEEG1{subji};
    
        

    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??


        % % %   % % Acquisition
        ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        %ids1 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        


        % % %   % % Acquisition only late trials
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) > 42;
        %ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) > 42;


        tfDCH1 = mean(EEG.power(ids1, :, : ,:), 'omitnan'); 
        tfDTF1 = squeeze(mean(tfDCH1, 2, 'omitnan'));
        tfDCH2 = mean(EEG.power(ids2, :, : ,:), 'omitnan'); 
        tfDTF2 = squeeze(mean(tfDCH2, 2, 'omitnan'));



        % % % % just to check number of trials
        %idsE1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) > 40;
        %idsL1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) <= 40;
        %xh1(subji) = sum(idsE1);         xh2(subji) = sum(idsL1); 
        %disp([num2str(sum(idsE1)) ' // ' num2str(sum(idsL1))])

        
        c1(subji, :, :) = tfDTF1; 
        c2(subji, :, :) = tfDTF2; 

                
    end

end


cd (paths.github)

%%

sub2exc = [];



c1B = c1(:, 3:8, 201:500); c2B = c2(:, 3:8, 201:500); 
%c1B = c1(:, 1:30, 201:500); c2B = c2(:, 1:30, 201:500); 
%c1B = c1(:, 1:54, :); c2B = c2(:, 1:54, :); 
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

%h = zeros(6, 300);
%h(clustinfo.PixelIdxList{4}) = 1; 

 





times = -1:.01:1.99; 
%times = -3:.01:3.99
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

%set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 ], 'yticklabels', {'1', '30'}, 'xlim', [-.5 2]);
set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power file2load '.png'], 'Resolution',150)
exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)







%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    %c1B = c1(:,1:30,301:480); 
    %c2B = c2(:,1:30,301:480); 
    c1B = c1(:,3:8,301:480); 
    c2B = c2(:,3:8,301:480); 
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :, :);
           c1B(subji, :, :) = c2B(subji, :, :);
           c2B(subji, :, :) = tmp; 
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





%% check only in theta band

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
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',150)
%% permutations 2D (line plot)

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = squeeze(mc1B); 
    c2B = squeeze(mc2B);
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



%% Correlate responses and power in cluster

clearvars -except ALLEEG ALLEEG1 paths clustinfo nL
close all
sub2exc = [];

for subji = 1:size(ALLEEG1, 1)
    EEG = ALLEEG1{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        %ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 201:500), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %cTR = squeeze(mean(powH(triali, :, 3:8, 321:390), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        nNan(subji,:) = sum(isnan(thPow));

        ratings2u = double(string(Ev2(ids, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        allNTR(subji,:) = length(thPow);



        allRho(subji, :) = corr(thPow, ratings2u, 'type', 'p');
        
        
        %figure()
        %scatter(thPow, ratings2u, 250, 'filled');
        %h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
        %C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
        %allSlopes(subji, :) = C(2);
        %allIntercepts(subji, :) = C(1);
        %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
        


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
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.85 .85] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


%% Extract data 4 LME

clearvars -except ALLEEG ALLEEG1 paths clustinfo nL
close all
sub2exc = [];

for subji = 1:size(ALLEEG1, 1)
    EEG = ALLEEG1{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        powH = EEG.power;
        
        for triali = 1:size(powH, 1)
            %cTR = squeeze(mean(powH(triali, :, 3:8, 201:500), 2));
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            cTR = squeeze(mean(powH(triali, :, 3:8, 301:350), 2));
            thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        ratings2u = double(string(Ev2(:, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        trialN = double(string(Ev2(:, 1)));
        exph = double(string(Ev2(:, 2)));
        trial_type = double(string(Ev2(:, 8)));

        trialN(ids2rem) = []; 
        trial_type(ids2rem) = []; 
        exph(ids2rem) = []; 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        subID = repelem(subji, length(ratings2u))';
        
        
        d4LME = [thPow ratings2u subID trialN exph trial_type];
        % % remove the testing phase for now
        d4LME(d4LME(:, 5) == 3,:) = []; 

        allData4LME{subji,:} = d4LME;
        
    end
end

%% 
d4LME = cat(1, allData4LME{:});

tbl2 = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), d4LME(:,5), d4LME(:,6), ...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'Phase', 'CS'});



%% fit model
clc


lme = fitlme(tbl2,'Ratings ~ theta_AMY + CS+ Phase+ trialN + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'Ratings ~ theta_AMY + CS+ Phase+ trialN + CS*theta_AMY + *theta_AMY + (1|subID)'); % random intercept model

lme
%lme.Coefficients














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

%EEG = ALLEEG{1}; 
tr = 20; 
ch = 2; 

% figure
% d2p	= squeeze(EEG.dsPower(tr, ch, : ,:));
% myCmap = colormap(brewermap([],'YlOrRd'));
% colormap(myCmap)
% contourf(1:300, 1:54, d2p, 40, 'linecolor', 'none'); colorbar

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
d2p	= squeeze(EEG.power(tr, ch, : ,:));
contourf(1:70, 1:54, d2p, 40, 'linecolor', 'none'); colorbar


%% plot example trial in one subject (ONLY 1 electrode)

%EEG = ALLEEG{1}; 
tr =26; 

figure
d2p	= squeeze(EEG.power(tr, 1 ,:,:));
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar


%% plot two different conditions (average in all amygdala channels and across trials)


Ev = [{EEG.event.type}]';
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
Ev2 = cat(1, Ev1{:})
ids = strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1') % CS+CS+ during acquisition
d2p1	= squeeze(mean(mean(EEG.power(ids, :, : ,:))));
ids = strcmp(Ev2(:, 6), '2')  & strcmp(Ev2(:, 2), '1') % CS+CS- during acquisition
d2p2	= squeeze(mean(mean(EEG.power(ids, :, : ,:))));

myCmap = colormap(brewermap([],'YlOrRd')); colormap(myCmap)

tiledlayout(2, 1,'TileSpacing','compact');
nexttile
contourf(1:701, 1:54, d2p1, 40, 'linecolor', 'none'); hold on; %colorbar
nexttile
contourf(1:701, 1:54, d2p2, 40, 'linecolor', 'none'); hold on; %colorbar












%% power analysis all channels / SETTINGS FOR PLOTTING IN BRAIN PLOTS 
clear, close all
paths = load_paths; 

c2u = 'C';
sROI = 'allChannels'; % don't change > for specific ROIs check next block 

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


for subji = 1:length(allsubs)
    subji
    sub = allsubs{subji}; 
    cd([paths.iEEG])
    load ([sub '_iEEG.mat']);

    
    Ev = [{EEG.event.type}]'; 
    Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
    clen = cellfun(@length, Ev1); 
    EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
    Ev2 = cat(1, Ev1{:});
   
    Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    ids = strcmp(Ev2(:, 10), c2u); 
    EEG.event = EEG.event(ids)

    %epoch data and markers
    EEG = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
    EEG = extract_power_EXT(EEG, 0.1); 
    EEG.power = EEG.power(:, :, :, 31:40);
    
    EEG = rmfield(EEG, 'data');
    EEG = normalize_EXT(EEG);

    ALLEEG{subji,:} = EEG; 


end




sROI = char(join(sROI, '_'));
filename = [paths.iEEGRes.power 'allS_' sROI '_' c2u];
save(filename, "ALLEEG");

cd (paths.github)


%% 