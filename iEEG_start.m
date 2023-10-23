%% 
%% 
 
clear, close all
paths = load_paths_EXT; 

c2u = 'C';

sROI = {'Amygdala'}; 
%sROI = {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'parahippocampal' 'entorhinal' };

%sROI = {'occipital' 'cuneus' 'lingual' 'pericalcarine' 'bankssts'}

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

        % % % remove renewal trials 
        %id2r3 = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2'); 
        %EEG.event = EEG.event(ids & id2r3);

        %epoch data and markers
        [EEG id2check] = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
        
        EEG = remove_elec_EXT(EEG, 50); %thres channels is 1/5 of 192 = 38

        

        if ~isempty(EEG.data)
            
            %EEG = extract_power_EXT(EEG, 0.01); 
            EEG = extract_theta_power_EXT(EEG); %HILBERT
            EEG = normalize_EXT(EEG);
            nChans(subji, :) = size(EEG.power, 2);
            ALLEEG{subji,:} = EEG; 

        end
    
    end


end

sROI = char(join(sROI, '_'));
mkdir(paths.results.power)
filename = [paths.results.power 'HILB_' sROI '_' c2u];
nSub = sum(cell2mat(cellfun(@(x) ~isempty(x), ALLEEG, 'un', 0)));
totalChans = sum(nChans);
save(filename, 'ALLEEG', 'nSub', 'nChans', 'totalChans', '-v7.3');

cd (paths.github)


%% % % % % % % % FIRST LOAD FILE - > START HERE
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 
%file2load = ['HILB_' 'Amygdala' '_C']; 

load ([paths.results.power file2load]); 



%% load traces

file2load = ['TR_' 'AMY' '_C']; 
load ([paths.results.traces file2load]); 


%% count channels 

clear nChans
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
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

%% PLOT grand average for each condition

clearvars -except ALLEEG ALLEEG1 paths  totalChans nChans nSub nL



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
        %ids3 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-


        % % %   % % Acquisition only late trials
        %ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & double(string(Ev2(:, 1))) > 42;
        %ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') & double(string(Ev2(:, 1))) > 42;


        tfDCH1 = mean(EEG.power(ids1, :, : ,:), 'omitnan'); 
        tfDTF1 = squeeze(mean(tfDCH1, 2, 'omitnan'));
        tfDCH2 = mean(EEG.power(ids2, :, : ,:), 'omitnan'); 
        tfDTF2 = squeeze(mean(tfDCH2, 2, 'omitnan'));
        %tfDCH3 = mean(EEG.power(ids3, :, : ,:), 'omitnan'); 
        %tfDTF3 = squeeze(mean(tfDCH3, 2, 'omitnan'));



        % % % % just to check number of trials
        %idsE1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) > 40;
        %idsL1 = ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1') & double(string(Ev2(:, 1))) <= 40;
        %xh1(subji) = sum(idsE1);         xh2(subji) = sum(idsL1); 
        %disp([num2str(sum(idsE1)) ' // ' num2str(sum(idsL1))])

        
        c1(subji, :, :) = tfDTF1; 
        c2(subji, :, :) = tfDTF2; 
        %c3(subji, :, :) = tfDTF3; 

                
    end

end


cd (paths.github)

%% ONLY THETA

sub2exc = [];



c1B = c1(:, 3:8, 201:500); c2B = c2(:, 3:8, 201:500); 
%c3B = c3(:, 3:8, 201:500); 
%c1B = c1(:, 1:30, 201:500); c2B = c2(:, 1:30, 201:500); 
%c1B = c1(:, 1:54, :); c2B = c2(:, 1:54, :); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 
%c3B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
%c3B(c3B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));
%d2p3	= squeeze(mean(c3B, 'omitnan'));


[h p ci ts] = ttest(c1B, c2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id)

% 

%h = zeros(6, 300);
%h(clustinfo.PixelIdxList{id}) = 1; 

 





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
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-3 3])
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


%% take in cluster only


for subji = 1:47

    c1BS = squeeze(c1B(subji,:,:));
    c2BS = squeeze(c2B(subji,:,:));
    c3BS = squeeze(c3B(subji,:,:));
    cSP(subji,:) = mean(c1BS(clustinfo.PixelIdxList{4}), 'all');
    cSMPP(subji,:) = mean(c2BS(clustinfo.PixelIdxList{4}), 'all');
    cSMPM(subji,:) = mean(c3BS(clustinfo.PixelIdxList{4}), 'all');

end


%% 
d2p = [cSP cSMPP cSMPM]

[~,~,stats] = anova1(d2p)
%[c,~,~,gnames] = multcompare(stats);

%%
clc

[h p ci ts] = ttest(cSP, cSMPP)
[h p ci ts] = ttest(cSP, cSMPM)

%%
d4ANOVA = [cSP; cSMPP ;cSMPM];
d4ANOVA(:,2) = [ones(1,47) ones(1,47)*2 ones(1,47)*3];
d4ANOVA(any(isnan(d4ANOVA), 2), :) = [];
d4ANOVA(:,3) = [1:32 1:32 1:32];

x = RMAOV1(d4ANOVA);

%% Correlate responses and power in cluster
load clustinfoAB
clearvars -except ALLEEG paths clustinfo nL pi2u
close all
sub2exc = []; %subj 23-35-38-39-44-45

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
        %ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 201:500), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        nNan(subji,:) = sum(isnan(thPow));

        ratings2u = double(string(Ev2(ids, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        allNTR(subji,:) = length(thPow);



        allRho(subji, :) = corr(thPow, ratings2u, 'type', 'k');
        
        
% % %         figure()
% % %         scatter(thPow, ratings2u, 250, 'filled');
% % %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% % %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% % %         allSlopes(subji, :) = C(2);
% % %         allIntercepts(subji, :) = C(1);
% % %         set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
        


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
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-1 1] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
set(gca, 'LineWidth', 3);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)







%% ALL FREQUENCIES

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

 





times = -.5:.01:1.79; 
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

set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'}, 'xlim', [-.5 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power file2load '.png'], 'Resolution',150)
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








%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    %c1B = c1(:,1:54,301:480); 
    %c2B = c2(:,1:54,301:480); 
    %c1B = c1(:,3:8,301:480); 
    %c2B = c2(:,3:8,301:480); 
    c1B = c1(:,3:8,301:400); % first second only
    c2B = c2(:,3:8,301:400); 
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
set(gca, 'xlim', [-.5 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',150)


%% THETA BAND - LINE PLOTS 2

sub2exc = [];

c1B = c1(:, 2001:5000); c2B = c2(:, 2001:5000); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
c1B(any(isnan(c1B), 2), :) = [];
c2B(any(isnan(c2B), 2), :) = [];

d2pm1	= squeeze(mean(c1B,'omitnan'));
d2pm2	= squeeze(mean(c2B,'omitnan'));
d2pstd1	= std(c1B);
d2pstd2	= std(c2B);
se1 = d2pstd1/sqrt(size(c1B, 1))
se2 = d2pstd2/sqrt(size(c2B, 1))

[h p ci ts] = ttest(c1B, c2B); 
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
times = (-1:.001:1.999) + .25;

colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 

%plot(times, d2p1); hold on; 
%plot(times, d2p2); hold on; 

xlabel('Time (s)')
ylabel('Theta Power')
plot (times, hb, 'Linewidth', 7)
%set(gca, 'xlim', [-.5 1.75])
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

load clustinfoAB
clearvars -except ALLEEG ALLEEG1 paths clustinfo nL pi2u
close all
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
            cTR = squeeze(mean(powH(triali, :, 3:8, 201:500), 2));
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            thPow(triali, :) = mean(cTR(pi2u), 'all');
            
            %cTR = squeeze(mean(powH(triali, :, 3:8, 301:470), 2));
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
        d4LME(d4LME(:, 5) == 3,:) = []; 

        allData4LME{subji,:} = d4LME;
        
    end
end

%% 
d4LME = cat(1, allData4LME{:});

%d2n = d4LME(:, 2); 
%mT = mean(d2n,'omitnan');
%stdT = std(d2n,[], 'omitnan');
%d4LME(:, 2) = bsxfun(@rdivide, bsxfun(@minus, d2n, mT), stdT);  
        


tbl2 = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), d4LME(:,5), d4LME(:,6), d4LME(:,7),...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'Phase', 'currCS', 'trial_type'});

tbl2.Ratings = ordinal(tbl2.Ratings);
%tbl2.currCS= categorical(tbl2.currCS);
tbl2.trial_type= categorical(tbl2.trial_type);
%tbl2.Phase= categorical(tbl2.Phase);


%%
VarDecompTbl = colldiag(d4LME)
% That can be passed along for visualization
colldiag_tableplot(VarDecompTbl);

%% fit model
clc


%lme = fitlme(tbl2,'Ratings ~ theta_AMY + trial_type + Phase + trialN + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'Ratings ~ theta_AMY + currCS+ Phase+ trialN + currCS*theta_AMY + Phase*theta_AMY + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'theta_AMY ~ Ratings + currCS + Phase+ trialN + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'Ratings ~ theta_AMY + currCS + Phase+ trialN + (1|subID)'); % random intercept model


lme = fitlme(tbl2,'theta_AMY ~ Ratings + currCS + Phase + Phase*currCS + (1|subID)'); % random intercept model

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