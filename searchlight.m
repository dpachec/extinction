%% LOAD ALLEEG, brain model and project data
%% 

paths = load_paths; 
file2load = 'allS_orbitofrontal_C'; 
%load ([paths.iEEGRes.power file2load]); 
clearvars -except ALLEEG paths file2load


c2u = file2load(end);


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        % % % % mean over trials for each electrode
        ids1 = strcmp(Ev2(:, 10), c2u) & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1');
        d2p1	= squeeze(mean(EEG.power(ids1, :, : ,:), 'omitnan'));
        ids2 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '3')  & strcmp(Ev2(:, 2), '1');
        d2p2	= squeeze(mean(EEG.power(ids2, :, : ,:), 'omitnan'));
    
    
        if ndims(d2p1) == 2 d2p12(1, :,:) = d2p1; d2p1 = d2p12; d2p22(1, :,:) = d2p2; d2p2 = d2p22; end %need this format later
        

        c1{subji, 1} = d2p1; 
        c1{subji, 2} = struct2cell(EEG.chanlocs)';
        c2{subji, 1} = d2p2; 
        c2{subji, 2} = struct2cell(EEG.chanlocs)';

        ch2u{subji,:} = struct2cell(EEG.chanlocs)';
        
    end

end


c1 = c1(~cellfun(@isempty, c1(:,1)), :);
c2 = c2(~cellfun(@isempty, c2(:,1)), :);
ch2u = ch2u(~cellfun(@isempty, ch2u(:,1)), :);

cd (paths.github)



%% load brain model 

%load surface
[vertices faces] = readObj('brain_LR.obj');

iDD = cellfun(@(x) cell2mat(x(:,3)),ch2u,'UniformOutput',false)


%% flip to use only 1 hemisphere

for subji = 1:length(chanNames_all)
    subjEle = iDD{subji};
    for ei = 1:size(subjEle, 1)
       if subjEle(ei,1) < 0
           einverted  = subjEle(ei,1)*-1;
           subjEle(ei,1) = einverted;
       end
    end
    iDD{subji} = subjEle;
end


%% find electrode idx

dist4elec = 15;

tic
clear eIds eIdist;
fprintf('\n'); fprintf('Subject:      '); %check this later
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    elec = iDD{subji};
    if ~isempty(elec)
        D = pdist2(elec,vertices)';
        t   = (D<dist4elec);
        for vi = 1:length(vertices)
           rowV = t(vi, :);
           rowVD = D(vi, :);
           idx = find(rowV);
           eIds{vi} = idx; % for each vertex the id of influencing elec
           s2check = rowVD(idx); 
           eIdist{vi} =  s2check; % for each vertex the distance of influencing elec
        end
        eIds = eIds'; 
        eIdist = eIdist';
    else
        eIds = [];
        eIdist = []; 
    end
    allEids{subji} = eIds;
    allEIdist{subji} = eIdist;
    
    clear eIds eIdist;
    
end

toc



%% count electrodes influencing each vertex
% subj2model is a subjects * vertices matrix containing average values to plot 
tic

clear subj2model
for subji = 1:length(iDD)
    for vxi = 1:length(vertices)
        if ~isempty(allEids{subji}{vxi}) 
            eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex

            d2proj1 = c1{subji, 1}(eId, :, :);
            d2proj2 = c2{subji, 1}(eId, :, :);
            d2proj = d2proj1-d2proj2; 

            % % % % project difference data
            d2projMeanChans = squeeze(mean(d2proj, 1)); % = freq by time
            d2projMeanChansFreq = squeeze(mean(d2projMeanChans(3:8,:), 1)); % only time
            d2projMeanChansFreqTime = squeeze(mean(d2projMeanChansFreq(11:20))); % all
            d2projF = d2projMeanChansFreqTime; 

            subj2model(subji, vxi) = d2projF;

            sumA(subji, vxi) = 1; % % % this is just to count subjects
        else
            subj2model(subji, vxi) = NaN;
        end
    end
end

toc

%% compute stats

minSubj = 6


t = squeeze(mean(subj2model, 'omitnan')); 


% % % % project statistics data 
%[h p ci ts] = ttest(subj2model); 
%t = ts.tstat; 


x = ~isnan(subj2model);
x1 = sum(x);
x2 = x1< minSubj; 
t(x2) = nan; 



toc



%% plot
clear gaPatch 

dat2plot        = 't'; %t or h
v2p             = 'l' %l or r
crange          = [-.05 .05];
%crange          = [-4 4];
%colorMap2use    = 'autumn';
colorMap2use    = 'jet';


gaPow2plot      = eval (dat2plot); %h or t or sumA
gaPatch = gaPow2plot; gaPatch(isnan(gaPow2plot)) = 0;  gaPatch(gaPatch ~= 0) = NaN;
tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch, gaPatch, gaPatch)';


figure(1);
vColor = ones(length(vertices), 3)-0.2;

pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines
pL2 = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_L,'FaceColor','interp'); hold on;
pL2.LineStyle = 'none';      % remove the lines

%pLZP = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_LZP,'FaceColor','interp');
%pLZP.LineStyle = 'none';      % remove the lines

caxis(crange)


%colorMap2use = 'hot';
colormap (colorMap2use); %winter works well (also try spring)
colorbar;
if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
else
	view(90,0)     % right orientation
end
%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .7 .3]) %sets the ambient/diffuse/specular strength of the objects.

set(gca, 'FontSize', 20)

axis equal off    % make the axes equal and invisible
%axis vis3d off
%axis manual
%p.FaceAlpha = 0.9;   % make the object semi-transparent
%pL.FaceColor = 'interp';    % set the face colors to be interpolated
%p.FaceColor = 'none';    % turn off the colors

%pR.LineStyle = 'none';      % remove the lines
x1= pL.FaceVertexCData;
finalColors = vals2colormap(x1, colorMap2use, crange);% NANs are converted to 1
finalColors (isnan(x1),:) = NaN;

%filename = ['_' dat2plot '_' num2str(f(1)) '-' num2str(f(end)) 'Hz.png'];
filename = 'elec_onset';
colormap(brewermap([],'*Spectral')); colorbar;
%brewermap_view({gca})
%export_fig(2, filename, '-r300', '-transparent');


%close all; 

% oldF = cd;
% cd '/Users/danielpacheco/Documents/UNITY/BX3_VR/Assets/MNI'
% csvwrite('colors.csv', finalColors);
% cd (oldF);
















%% plot electrodes frontal view

cd (paths.elec_images)

nSubj = length(iDD);

radius = 1.8; 
cols = brewermap(length(iDD), 'Dark2'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
close all
lighting none

figure(1); set(gcf, 'Position', [100 100 1500 1500])
tiledlayout(6, 8,'TileSpacing','none', 'Padding','none');
for subji = 1:nSubj
    nexttile
    vColor = ones(length(vertices), 3)-0.2;
    pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
    pL.LineStyle = 'none';      % remove the lines
    alpha 0.4
    view(180,0)
    %campos([1,1,1])
    %l = light('Position',[0.9 1 -1],'Style','infinite')

    count = 1;

    
    clear subjEle
    subjEle = iDD{subji};

    %use this to flip electrodes from one hemisphere 
%     for ei = 1:length(subjEle)
%        if subjEle(ei,1) > 0 %if they are on left (or right) hemisphere change
%            einverted  = subjEle(ei,1)*-1;
%            subjEle(ei,1) = einverted;
%        end
%     end
    
    for ei = 1:size(subjEle, 1)
        [x, y, z] = sphere;
        x = x*radius; y = y*radius; z = z*radius; 
        mesh(x + subjEle(ei, 1), y +  subjEle(ei, 2) , z + subjEle(ei, 3), 'edgecolor', 'none', 'facecolor', cols(count,:));
    end
    count = count+1; 
    set(gca,'XTick',[],'YTick',[])
    set(gca, 'Visible', 'off')
    axis equal off    % make the axes equal and invisible
end


        
    filename = [ 'all_elec_.png'];
    exportgraphics(gcf, filename,'Resolution',300)
    %close all




 




close all
disp('done');





%% plot electrodes left view

cd (paths.elec_images)


radius = 3; 
cols = brewermap(length(iDD), 'Dark2'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
close all
lighting none

figure(1); set(gcf, 'Position', [100 100 1200 1200])
tiledlayout(8, 6,'TileSpacing','none', 'Padding','none');
for subji = 1:nSubj
    nexttile
    vColor = ones(length(vertices), 3)-0.2;
    pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
    pL.LineStyle = 'none';      % remove the lines
    alpha 0.3
    view(270,0)
    x = [-110 80]; 
    y = [-55 85]; 
    [x y] = meshgrid(x, y); %
    z = zeros(size(x, 1)); % Generate z data
    h = mesh(z,x,y); % Plot the surface
    set(h,'edgecolor','none','facecolor',[1 1 1 ])
    
    count = 1;

    
    clear subjEle
    subjEle = iDD{subji};

    %use this to flip electrodes from one hemisphere 
%     for ei = 1:length(subjEle)
%        if subjEle(ei,1) > 0 %if they are on left (or right) hemisphere change
%            einverted  = subjEle(ei,1)*-1;
%            subjEle(ei,1) = einverted;
%        end
%     end
    
    for ei = 1:size(subjEle, 1)
        [x, y, z] = sphere;
        x = x*radius; y = y*radius; z = z*radius; 
        mesh(x + subjEle(ei, 1), y +  subjEle(ei, 2) , z + subjEle(ei, 3), 'edgecolor', 'none', 'facecolor', cols(count,:));
    end
    count = count+1; 
    set(gca,'XTick',[],'YTick',[])
    set(gca, 'Visible', 'off')
    axis equal off    % make the axes equal and invisible
end


        
    filename = [ 'all_elec_left.png'];
    exportgraphics(gcf, filename,'Resolution',300)

%close all
disp('done');

%% plot electrodes right view

cd (paths.elec_images)


radius = 3; 
cols = brewermap(length(iDD), 'Dark2'); % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
close all
lighting none

figure(1); set(gcf, 'Position', [100 100 1200 1200])
tiledlayout(8, 6,'TileSpacing','none', 'Padding','none');
for subji = 1:nSubj
    nexttile
    vColor = ones(length(vertices), 3)-0.2;
    pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
    pL.LineStyle = 'none';      % remove the lines
    alpha 0.3
    view(90,0)
    x = [-110 80]; 
    y = [-55 85]; 
    [x y] = meshgrid(x, y); %
    z = zeros(size(x, 1)); % Generate z data
    h = mesh(z,x,y); % Plot the surface
    set(h,'edgecolor','none','facecolor',[1 1 1 ])
    
    count = 1;

    
    clear subjEle
    subjEle = iDD{subji};

    %use this to flip electrodes from one hemisphere 
%     for ei = 1:length(subjEle)
%        if subjEle(ei,1) > 0 %if they are on left (or right) hemisphere change
%            einverted  = subjEle(ei,1)*-1;
%            subjEle(ei,1) = einverted;
%        end
%     end
    
    for ei = 1:size(subjEle, 1)
        [x, y, z] = sphere;
        x = x*radius; y = y*radius; z = z*radius; 
        mesh(x + subjEle(ei, 1), y +  subjEle(ei, 2) , z + subjEle(ei, 3), 'edgecolor', 'none', 'facecolor', cols(count,:));
    end
    count = count+1; 
    set(gca,'XTick',[],'YTick',[])
    set(gca, 'Visible', 'off')
    axis equal off    % make the axes equal and invisible
end


        
    filename = [ 'all_elec_right.png'];
    exportgraphics(gcf, filename,'Resolution',300)

%close all
disp('done');


%% find electrode idx

dist4elec = 12.5;

tic
clear eIds eIdist;
fprintf('\n'); fprintf('Subject:      '); %check this later
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    elec = iDD{subji};
    if ~isempty(elec)
        D = pdist2(elec,vertices)';
        t   = (D<dist4elec);
        for vi = 1:length(vertices)
           rowV = t(vi, :);
           rowVD = D(vi, :);
           idx = find(rowV);
           eIds{vi} = idx; % for each vertex the id of influencing elec
           s2check = rowVD(idx); 
           eIdist{vi} =  s2check;
        end
        eIds = eIds'; % for each vertex the id of influencing elec
        eIdist = eIdist';
    else
        eIds = [];
        eIdist = []; 
    end
    allEids{subji} = eIds;
    allEIdist{subji} = eIdist;
    
    clear eIds eIdist;
    
end

toc


%% count electrodes influencing each vertex

tic
minSubjects = 1;

fprintf('\n'); fprintf('Subject:      '); %check this later
clear percContr2plot
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    for vxi = 1:length(vertices)
        if ~isempty(allEids{subji}{vxi}) 
            eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex
            percContr2plot(subji, vxi) = 1;
        else
            percContr2plot(subji, vxi) = NaN;
        end

    end
        
end


idx2c = percContr2plot; %idx2c = for each vertex, presence or absence of subjects
sumA = sum(idx2c, 'omitnan'); % sumA = number of subjects affecting each vertex
idx2c2 = find (sumA < minSubjects); % ids of vtx with less than minSubj
sumA(:, idx2c2)= 0;
%percContr2plot(:, idx2c2) = NaN; 
%[h p ci t] = ttest (pow2Plot_AC, pow2Plot_PC, 'Alpha', 0.05); % grand average 
%t = mean(percContr2plot, 'omitnan'); %t.tstat;
t = sumA; 

toc

%% plot

dat2plot        = 't'; %t or h
v2p             = 'r' %l or r
crange          = [1 12];%[-.3 .5];
%crange         = [-10 10];
cmp2use         = 'YlOrRd'%'YlOrRd';

gaPow2plot      = eval (dat2plot); %h or t or sumA
gaPatch = gaPow2plot; gaPatch(isnan(gaPow2plot)) = 0;  gaPatch(gaPatch ~= 0) = NaN;
% subj = 13;
% tmean_L = pow2Plot(subj,:);
% tmean_LZP = cat (1, zero2Plot(subj,:),  zero2Plot(subj,:),  zero2Plot(subj,:));

%tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch+.5, gaPatch+.5, gaPatch+.5)';
tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch, gaPatch, gaPatch)';


figure(2);
%cindex = cat(1, tmean, tmean, tmean)';
%cindex = ones(length(g.vertices), 3);
%cindex(t1) = zeros(1,3);
pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_L,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines


%lighting gouraud; %p.SpecularStrength = 0.4;%camlight

pLZP = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_LZP,'FaceColor','interp');
pLZP.LineStyle = 'none';      % remove the lines
%pR = patch('Faces',gR.faces,'Vertices',gR.vertices,'FaceVertexCData',tmean_R','FaceColor','interp')

caxis(crange)


colorMap2use = 'hot';
colormap (colorMap2use); %winter works well (also try spring)
colorbar;
if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
else
	view(90,0)     % right orientation
end
%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .7 .3]) %sets the ambient/diffuse/specular strength of the objects.

set(gca, 'FontSize', 20)

axis equal off    % make the axes equal and invisible
%axis vis3d off
%axis manual
%p.FaceAlpha = 0.9;   % make the object semi-transparent
%pL.FaceColor = 'interp';    % set the face colors to be interpolated
%p.FaceColor = 'none';    % turn off the colors

%pR.LineStyle = 'none';      % remove the lines
x1= pL.FaceVertexCData;
finalColors = vals2colormap(x1, colorMap2use, crange);% NANs are converted to 1
finalColors (isnan(x1),:) = NaN;

%filename = ['_' dat2plot '_' num2str(f(1)) '-' num2str(f(end)) 'Hz.png'];
colormap(brewermap([],cmp2use)); colorbar; %PuBuGn
%brewermap_view({gca})

filename = [v2p '_t_contribution.jpg']
export_fig(2, filename, '-r300', '-transparent');


%close all; 

% oldF = cd;
% cd '/Users/danielpacheco/Documents/UNITY/BX3_VR/Assets/MNI'
% csvwrite('colors.csv', finalColors);
% cd (oldF);



%% find electrode idx

dist4elec = 20;

tic
clear eIds eIdist;
fprintf('\n'); fprintf('Subject:      '); %check this later
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    elec = iDD{subji};
    if ~isempty(elec)
        D = pdist2(elec,vertices)';
        t   = (D<dist4elec);
        for vi = 1:length(vertices)
           rowV = t(vi, :);
           rowVD = D(vi, :);
           idx = find(rowV);
           eIds{vi} = idx; % for each vertex the id of influencing elec
           s2check = rowVD(idx); 
           eIdist{vi} =  s2check;
        end
        eIds = eIds'; % for each vertex the id of influencing elec
        eIdist = eIdist';
    else
        eIds = [];
        eIdist = []; 
    end
    allEids{subji} = eIds;
    allEIdist{subji} = eIdist;
    
    clear eIds eIdist;
    
end

toc
%% assign pow4percContr values to vertices 
tic
minSubjects = 5;
f = [3:8];
percContrPow = mean(pow4percContr(:,:, f), 3);

perc2Plot = zeros (length(allEids), length(vertices));
fprintf('\n'); fprintf('Subject:      '); %check this later
clear percContr2plot
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    perCvalues = []; 
    for vxi = 1:length(vertices)
        if ~isempty(allEids{subji}) & ~isempty(percContrPow(:, subji))
            eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex
            
            allM  = percContrPow(eId,subji);
            perCvalues = mean(allM);
            
        
            percContr2plot(subji, vxi) = perCvalues;
            
        else
            percContr2plot(subji, vxi) = NaN;
        end

    end
        
end


idx2c = ~isnan(percContr2plot); %idx2c = for each vertex, presence or absence of subjects
sumA = sum(idx2c, 'omitnan'); % sumA = number of subjects affecting each vertex
idx2c2 = find (sumA < minSubjects); % ids of vtx with less than minSubj
sumA(:, idx2c2)= NaN;
percContr2plot(:, idx2c2) = NaN; 
[h p ci t] = ttest (percContr2plot); % grand average 
%t = mean(percContr2plot, 'omitnan'); %t.tstat;
t = t.tstat;

toc





%% plot

dat2plot        = 'h'; %t or h
v2p             = 'l' %l or r
crange          = [-1 4];
%colorMap2use    = 'autumn';
colorMap2use    = 'jet';


gaPow2plot      = eval (dat2plot); %h or t or sumA
gaPatch = gaPow2plot; gaPatch(isnan(gaPow2plot)) = 0;  gaPatch(gaPatch ~= 0) = NaN;
% subj = 13;
% tmean_L = pow2Plot(subj,:);
% tmean_LZP = cat (1, zero2Plot(subj,:),  zero2Plot(subj,:),  zero2Plot(subj,:));

%tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch+.5, gaPatch+.5, gaPatch+.5)';
tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch, gaPatch, gaPatch)';


figure(2);
%cindex = cat(1, tmean, tmean, tmean)';
%cindex = ones(length(g.vertices), 3);
%cindex(t1) = zeros(1,3);
pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_L,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines


%lighting gouraud; %p.SpecularStrength = 0.4;%camlight

pLZP = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_LZP,'FaceColor','interp');
pLZP.LineStyle = 'none';      % remove the lines
%pR = patch('Faces',gR.faces,'Vertices',gR.vertices,'FaceVertexCData',tmean_R','FaceColor','interp')

caxis(crange)


%colorMap2use = 'hot';
colormap (colorMap2use); %winter works well (also try spring)
colorbar;
if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
else
	view(90,0)     % right orientation
end
%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .7 .3]) %sets the ambient/diffuse/specular strength of the objects.

set(gca, 'FontSize', 20)

axis equal off    % make the axes equal and invisible
%axis vis3d off
%axis manual
%p.FaceAlpha = 0.9;   % make the object semi-transparent
%pL.FaceColor = 'interp';    % set the face colors to be interpolated
%p.FaceColor = 'none';    % turn off the colors

%pR.LineStyle = 'none';      % remove the lines
x1= pL.FaceVertexCData;
finalColors = vals2colormap(x1, colorMap2use, crange);% NANs are converted to 1
finalColors (isnan(x1),:) = NaN;

filename = ['_' dat2plot '_' num2str(f(1)) '-' num2str(f(end)) 'Hz.png'];
%filename = 'elec_onset';
colormap(brewermap([],'*Spectral')); colorbar;
%brewermap_view({gca})
export_fig(2, filename, '-r300', '-transparent');


%close all; 

% oldF = cd;
% cd '/Users/danielpacheco/Documents/UNITY/BX3_VR/Assets/MNI'
% csvwrite('colors.csv', finalColors);
% cd (oldF);




%% assign elecOnset values to vertices 
tic
minSubjects = 2;

perc2Plot = zeros (length(allEids), length(vertices));
fprintf('\n'); fprintf('Subject:      '); %check this later
clear percContr2plot
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    perCvalues = []; 
    for vxi = 1:length(vertices)
        if ~isempty(allEids{subji}) & ~isempty(elecOnset(:, subji))
            eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex
            
            allM  = elecOnset(eId,subji);
            perCvalues = mean(allM);
            
        
            percContr2plot(subji, vxi) = perCvalues;
            
        else
            percContr2plot(subji, vxi) = NaN;
        end

    end
        
end


idx2c = ~isnan(percContr2plot); %idx2c = for each vertex, presence or absence of subjects
sumA = sum(idx2c, 'omitnan'); % sumA = number of subjects affecting each vertex
idx2c2 = find (sumA < minSubjects); % ids of vtx with less than minSubj
sumA(:, idx2c2)= NaN;
percContr2plot(:, idx2c2) = NaN; 
%[h p ci t] = ttest (percContr2plot); % grand average 
t = mean(percContr2plot, 'omitnan'); %t.tstat;
%t = t.tstat;

toc




%% electrodes 2 use from excel and Channames based on those with (power)values
%note that electrodes are sorted differently in allElecM and chanNames_all

load chanNames_all;
clear matchEl allElecM elFromExc;
for subji = 1:13
    subji
    elfromRSA  = chanNames_all{subji};
    %elwithPowV = powerAllElec{subji, 2};
    elwithPowV = chanNames_all{subji};
    [elwithPowandRSA idmE{subji}] = intersect(elfromRSA, elwithPowV ); 
    elFromExcl  = allElec{subji, 1}(:,1);
    [matchEl{subji} idmE{subji}] = intersect(elFromExcl, elwithPowandRSA ); 
    allElecM{subji} = allElec{subji}(idmE{subji}, :);
   
end
matchEl = matchEl';
idmE = idmE';
allElecM = allElecM';



%% get the MNI coordinates from excel file
tic 
clear i index allReg iD iDD;
for subji = 1:13
   
        clear c x;
        c = allElecM{subji};
        d = c(:,3);
        iD = string(d);
        if size(iD, 1) > 0
            for j = 1:length(iD) 
               x(j,:) = str2num(iD(j)); 
            end
        else
            disp ('no x');
            x =[];
        end
        iDD{subji} = x;  
    
end
iDD = iDD';


toc

%% plot coverage

v2p = 'l' %l or r

cols = jet(13);


figure(2);
%cindex = cat(1, tmean, tmean, tmean)';
vColor = ones(length(vertices), 3)-0.15;
%cindex(t1) = zeros(1,3);
pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines
alpha 0.2

for subji = 1:length(iDD)
    subji
    subjEle = iDD{subji};
    for ei = 1:length(subjEle)
        [x, y, z] = sphere;
        mesh(x + subjEle(ei, 1), y +  subjEle(ei, 2) , z + subjEle(ei, 3), 'edgecolor', cols(subji,:));
        
    end
    
end



if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
else
	view(90,0)     % right orientation
end
%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .3 .3]) %sets the ambient/diffuse/specular strength of the objects.


axis equal off    % make the axes equal and invisible
%axis vis3d off
%axis manual
%p.FaceAlpha = 0.9;   % make the object semi-transparent
%pL.FaceColor = 'interp';    % set the face colors to be interpolated
%p.FaceColor = 'none';    % turn off the colors

%pR.LineStyle = 'none';      % remove the lines

export_fig(2, '_electrode_coverage_L.png', '-r300', '-transparent');
close all;

%% assign power values to vertices ALL TRIALS
tic
f           = [39:54];
minSubjects = 1;

pow2Plot = zeros (length(allEids), length(vertices));
fprintf('\n'); fprintf('Subject:      '); %check this later
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    powValues = [];
    for vxi = 1:length(vertices)

        eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex
        
        all  = powerAllElec{subji}(:, eId, :);
        %isempty(all)

        if ~isempty(all)
            powValuesAcT = mean (all, 1);
            powValuesAcF = mean(powValuesAcT(:,:,f), 3);
            powValues = squeeze(mean(powValuesAcF));
            pow2Plot(subji, vxi) = powValues;
        end

        
    end
        
end


idx2c = ~isnan(pow2Plot); %idx2c = for each vertex, presence or absence of subjects
sumA = sum(idx2c, 'omitnan'); % sumA = number of subjects affecting each vertex
idx2c2 = find (sumA < minSubjects); % ids of vtx with less than minSubj
pow2Plot_C = pow2Plot; 
pow2Plot_C(:, idx2c2) = NaN; 
sumA(:, idx2c2) = NaN;
t = mean(pow2Plot_C, 'omitnan'); %mean across subjects

toc


%% plot

dat2plot        = 't'; %t or h
v2p             = 'l' %l or r
crange          = [0 .15 ];
colorMap2use    = 'autumn';

gaPow2plot      = eval (dat2plot); %h or t or sumA
gaPatch = gaPow2plot; gaPatch(isnan(gaPow2plot)) = 0;  gaPatch(gaPatch ~= 0) = NaN;
% subj = 13;
% tmean_L = pow2Plot(subj,:);
% tmean_LZP = cat (1, zero2Plot(subj,:),  zero2Plot(subj,:),  zero2Plot(subj,:));

%tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch+.5, gaPatch+.5, gaPatch+.5)';
tmean_L = gaPow2plot';  tmean_LZP = cat (1, gaPatch, gaPatch, gaPatch)';


figure(2);
%cindex = cat(1, tmean, tmean, tmean)';
%cindex = ones(length(g.vertices), 3);
%cindex(t1) = zeros(1,3);
pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_L,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines


%lighting gouraud; %p.SpecularStrength = 0.4;%camlight

pLZP = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',tmean_LZP,'FaceColor','interp');
pLZP.LineStyle = 'none';      % remove the lines
%pR = patch('Faces',gR.faces,'Vertices',gR.vertices,'FaceVertexCData',tmean_R','FaceColor','interp')

caxis(crange)


%colorMap2use = 'hot';
colormap (colorMap2use); %winter works well (also try spring)
colorbar;
if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
else
	view(90,0)     % right orientation
end
%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .7 .3]) %sets the ambient/diffuse/specular strength of the objects.


axis equal off    % make the axes equal and invisible
%axis vis3d off
%axis manual
%p.FaceAlpha = 0.9;   % make the object semi-transparent
%pL.FaceColor = 'interp';    % set the face colors to be interpolated
%p.FaceColor = 'none';    % turn off the colors

%pR.LineStyle = 'none';      % remove the lines
x1= pL.FaceVertexCData;
finalColors = vals2colormap(x1, colorMap2use, crange);% NANs are converted to 1
finalColors (isnan(x1),:) = NaN;

filename = ['_' dat2plot '_' num2str(f(1)) '-' num2str(f(end)) 'Hz.png'];
%filename = ['image.jpg']
export_fig(2, filename, '-r300');
close all; 

% oldF = cd;
% cd '/Users/danielpacheco/Documents/UNITY/BX3_VR/Assets/MNI'
% csvwrite('colors.csv', finalColors);
% cd (oldF);




%% assign power values to vertices ACTIVE VS PASSIVE
tic
f           = [30:54];
minSubjects = 3;

pow2Plot = zeros (length(allEids), length(vertices));
fprintf('\n'); fprintf('Subject:      '); %check this later
for subji = 1:length(iDD)
    if (subji < 10) fprintf('\b\b\b'); elseif (subji < 100) fprintf('\b\b\b\b'); end
    fprintf('%d %s', subji, ' '); 
    powValues = [];
    for vxi = 1:length(vertices)
        if ~isempty(allEids{subji}) & ~isempty(powerAllElec{subji, 1})
            eId = cell2mat(allEids{subji}(vxi)); %elec id for each vertex
            %eIdist = cell2mat(allEids{subji}(vxi));
            
            all  = powerAllElec{subji, 1}(eId)';
            allM = cell2mat(powerAllElec{subji, 1}(eId));
            if length(all) > 1
                %disp ('reshape');
                powValues = reshape (allM, [length(all), 100, 6]); % change to 2 or 6 if sme
                %disp (size (powValues));
                
            else
                %disp ('no reshape');
                powValues = allM;
                %disp (size (powValues));
            end
            %powVal2check{vxi} = powValues;
        end 
        if ~isempty(powValues)
            %average from all influencing electrodes
            %size(powValues) %powValues has the shape (nEl_influencing, 100, 6)
            if ndims(powValues) == 2 % >> no reshape
                %disp ('only 1 electrode for this vertex');
                %only one influencing electrode
                powVBand_A = mean (powValues(f,1), 1); % in a specific band
                powVBand_P = mean (powValues(f,2), 1); % in a specific band
                
            else
                %disp ('more than 1 electrode for this vertex');
                meanPV = mean(powValues, 'omitnan'); %mean across influencing electrodes
                powMeanE = squeeze(meanPV); 
                
                powVBand_A = mean (powMeanE(f,1), 1); % in a specific band
                powVBand_P = mean (powMeanE(f,2), 1); % 
                
            end
            
        else
            powVBand_A = NaN;
            powVBand_P = NaN;
        
        end
        
        pow2Plot_A(subji, vxi) = powVBand_A;
        pow2Plot_P(subji, vxi) = powVBand_P;

    end
        
end


idx2c = ~isnan(pow2Plot_A); %idx2c = for each vertex, presence or absence of subjects
sumA = sum(idx2c, 'omitnan'); % sumA = number of subjects affecting each vertex
idx2c2 = find (sumA < minSubjects); % ids of vtx with less than minSubj
pow2Plot_AC = pow2Plot_A; pow2Plot_PC = pow2Plot_P;
pow2Plot_AC(:, idx2c2) = NaN; pow2Plot_PC(:, idx2c2) = NaN;
sumA(:, idx2c2) = NaN;
[h p ci t] = ttest (pow2Plot_AC, pow2Plot_PC, 'Alpha', 0.05); % grand average 
%t = mean(pow2Plot_PC, 'omitnan'); %t.tstat;
t = mean(pow2Plot_PC, 'omitnan'); %t.tstat;

toc




%% take hippocampal as a validity check

H = [10, 1, 7, 1, 11, 20, 15, 10, 10, NaN, 4, 8, NaN];
%H = [12, 12, 17, 12, 21, 22, 5, 7, 2, 4, 5, 7, 9];

clear powE powE_A powE_P;
for subji = 1:13
    subji
    if ~isnan(H(subji))
        powE = cell2mat(powerAllElec{subji,1}(H(subji)));
        powE_A(:, subji) = powE(:,1);
        powE_P(:, subji) = powE(:,2);
    end
    
end

mA = mean(powE_A, 2);
mP = mean(powE_P, 2);

figure()
plot (1:100, mA, 'r'); hold on;
plot (1:100, mP, 'b');

[h p ci t] = ttest(powE_A', powE_P');







%% load gii files
cd '/Users/danielpacheco/Documents/BrainX3vis/bids/sub-mni/anat'
g = gifti('sujet01_Lwhite.surf.gii')


%% load gifty gii
tic     
%clearvars -except ALLEEG
clearvars 

sublist = dir('*.surf.gii');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    g = gifti(sublist{subji});
    filename = [sublist{subji} 'obj'];
    saveas(g, filename, 'obj');
    %s_all{subji} = export(g);
end
s_all = s_all';
clearvars -except s_all sublist
toc
%% convert percContr to pow_all_elec sme


for subji = 1:13
    
    perc_all_elec = percContr
    
    
    
end 





%%
g = s_all{1};


%%

EEG = ALLEEG{1};
data = EEG.data(48,:);

data1 = downsample(data, 10);
mean_data = mean(data1, 'omitnan'); std_data = std(data1, [], 2, 'omitnan');
data2 = bsxfun(@rdivide, bsxfun(@minus, data, mean_data), std_data);
myM = hot(256);
data3 = vals2colormap(data2,  'hot', [-3 3]);

csvwrite('test1.csv', data3);



%% plot china data
%firs load obj 
[vertices faces] = readObj('brain_LR.obj');

%%
sublist = dir('*_frv.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load (sublist{subji});
    clear myCoord;
    subjEleAll{subji} = elec_mni_frv.chanpos;
end

v2p = 'l' %l or r
cols = jet(6);

figure();
vColor = ones(length(vertices), 3)-0.15;
pL = patch('Faces',faces,'Vertices',vertices,'FaceVertexCData', vColor,'FaceColor','interp'); hold on;
pL.LineStyle = 'none';      % remove the lines
alpha 0.2

for subji = 1:6
    subjEle = subjEleAll{subji};
    for ei = 1:length(subjEle)
        [x, y, z] = sphere;
        mesh(x + subjEle(ei, 1), y +  subjEle(ei, 2) , z + subjEle(ei, 3), 'edgecolor', cols(subji,:));
        
    end
    
end



if (strcmp (v2p, 'l'))
    view(270,0)     % left orientation
    view(180,0)    
else
	view(90,0)     % right orientation
    view(180,0)    
end

%set(gca, 'CameraViewAngle', 3); %for orthographic positionnig. 
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .3 .3]) %sets the ambient/diffuse/specular strength of the objects.


axis equal off    % make the axes equal and invisible





%%
%% LOAD ELECTRODE DATA, SURFACE AND POWER RESULTS OBJ
%(run within all_electrodes folder)
clearvars %-except chanNames_all
tic


%first load chinese data
cd ('C:\Users\Neuropsychology\Desktop\_WM\electrode_MNI\china')

sublist = dir('*_frv.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load (sublist{subji});
    clear myCoord;
    subjEleAll{subji} = elec_mni_frv.chanpos;
end


cd .. 


load ('chanNames_all')

ch2u{1} = chanNames_all{1};
ch2u{2} = chanNames_all{2};
ch2u{3} = chanNames_all{4};
ch2u{4} = chanNames_all{6};
ch2u{5} = chanNames_all{8};
ch2u{6} = chanNames_all{10};
ch2u{7} = chanNames_all{15};
ch2u{8} = chanNames_all{19};
ch2u{9} = chanNames_all{21};
ch2u{10} = chanNames_all{23};
ch2u{11} = chanNames_all{25};
ch2u{12} = chanNames_all{27};

clear chanNames_all
ch2u = ch2u'; 

%load surface
%[vertices faces] = readObj('low_res_L&R90.obj');
[vertices faces] = readObj('brain_LR.obj');

%load ('pow4percContr_R');
%load ('power_all_elec_sme_1-54Hz');
%load ('power_all_elec_1_100');
%load ('percContr_SI_norm');
%load ('percContr_SIA_norm');

%load ('phaseAllElec'); powerAllElec = phaseAllElec;

S = load ('all_elec'); %.mat file from pipEEG
S = struct2cell(S);
allElecReal = S{1};

%%sort

for subji = 1:12
   subji
   clear array2 reference x
   reference = ch2u{subji};
   array2 = allElecReal{subji, 2};
   [~, pos] = ismember(reference,array2); 
   pos(pos == 0) = [];
   allElecReal1{subji, 2} = array2(pos);
   x = allElecReal{subji, 1}(pos,:);
   allElecReal1{subji, 1} = x;
end

iDD = allElecReal1(:,1) ;


iDD(13:18) = subjEleAll;


disp('sorted');


%%

nL = cellfun(@(x) size(x,1), chanNames_all, 'un', 0);
c2w = cat(1, chanNames_all{:})
T = cell2table(c2w)
writetable(T,'myDataFile.csv')








%% count channels in each region 
clear
tic 
paths = load_paths
cd (paths.elec)


load allSChans
ch2u = allSChans'; 

for subji = 1:length(allSChans)
    
    elecS = allSChans{subji}    ;
    totalNofElec(subji,:) = length(elecS);
    chansLab = elecS(:, 2);
    selChans = contains(chansLab, 'Amygdala');
    sCH_Amygdala{subji,:} = sum(selChans);
    selChans = contains(chansLab, 'Hippocampus');
    sCH_Hippocampus{subji,:} = sum(selChans);

    

end

%%

sum(totalNofElec)
sum(cell2mat(sCH_Hippocampus))
sum(cell2mat(sCH_Amygdala))



cell2mat(sCH_Amygdala) == 0
















%%













