%% % % load edf and trial info file

eeglab
clear, close all

load /Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/iEEG/preproc/trialinfo/c_sub01_trlinfo.mat
%load D:/extinction/data/preproc/trialinfo/c_sub02_trlinfo.mat
%load D:/extinction/data/preproc/trialinfo/c_sub03_trlinfo.mat

EEG = pop_biosig('/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/raw_data/c_sub01/ieeg/DBX-Learning test.edf', 'importevent', 'off'); 
%EEG = pop_biosig('D:\extinction\raw_data\c_sub02\ieeg\LSY-Learning test.edf', 'importevent', 'off'); 
%EEG = pop_biosig('D:\extinction\raw_data\c_sub03\ieeg\JJL-Learning test.edf', 'importevent', 'off'); 

%% Extract triggers from TTL channel
eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X > 100000';

EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on');

chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 100, 'spacing', 10000000, 'events', EEG.event);


%%
factor = ( 10000 / EEG.srate); 

x = [EEG.event.latency]';
[d id2Break] = max(abs(diff(x)));
plot(abs(diff(x)))
lat2breakTTL = EEG.event(id2Break+1).latency; 


lat_prev = trlinfo(:, 10:14); 
sCod = string(trlinfo(:, 1:9));
sCod1 = repelem(sCod, 5, 1);
sCodF = join(sCod1, '_');

latencies = reshape (lat_prev', 1, [])'; 
%latencies = latencies(~isnan(latencies));


x = latencies;
[d id2Break] = max(abs(diff(x)));
%plot(diff(x))
%lat2breakLog = 1444000; %x(id) ; 
lat2breakLog = x(id2Break+2)/factor;
newBreak = lat2breakTTL - lat2breakLog


%diff_13 = EEG.event(id+1).latency;

%% take second 

factor = ( 10000 / EEG.srate); 
lat = EEG.event(1).latency; 
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
diff_12 =  latencies(2) / factor - lat;

for i = 1:id2Break
    EEG.event(i) = struct('latency', latencies(i)/factor -diff_12 , 'type', sCodF(i,:), 'urevent', 0);
end
for i = id2Break+1:length(latencies)
     EEG.event(i) = struct('latency', latencies(i)/factor + newBreak , 'type',sCodF(i,:), 'urevent', 0);
end


eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);




%% ONE AT THE TIME 

%eeglab
clear, close all
paths = load_paths; 
currentPath = pwd; 

sub = 'c_sub29'

cd(paths.trlinfo)
log_list = dir('*mat'); log_list = {log_list.name};

cd(paths.raw_data)
edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
edf_list = {edf_list.name}';

cd(paths.trlinfo)
load ([sub '_trlinfo.mat']);

cd([paths.raw_data sub '/ieeg/'])
edf_list = dir('*edf'); edf_list = {edf_list.name};
if length(edf_list) < 2    
    EEG = pop_biosig(edf_list{1}, 'importevent', 'off'); 
else
    EEG = pop_biosig(edf_list{1}, 'importevent', 'off'); 
    EEG2 = pop_biosig(edf_list{2}, 'importevent', 'off'); EEG.data = [EEG.data EEG2.data];
end        

clearvars -except EEG subji paths trlinfo path_log currentPath sub
eventChannel = 'POL DC12';
if strcmp(sub, 'c_sub09') eventChannel = 'POL DC09'; end

TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X > 100000';
EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on');

factor = ( 10000 / EEG.srate); 
x = [EEG.event.latency]';
[d id2Break] = max(abs(diff(x)));

if strcmp(sub, 'c_sub04') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-2);end
if strcmp(sub, 'c_sub06') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-1);end
if strcmp(sub, 'c_sub24') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-4);end
if strcmp(sub, 'c_sub25') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-8);end

lat2breakTTL = EEG.event(id2Break+1).latency;
if strcmp(sub, 'c_sub09') lat2breakTTL = EEG.event(id2Break+10).latency; end
if strcmp(sub, 'c_sub12') lat2breakTTL = EEG.event(id2Break+63).latency; end
if strcmp(sub, 'c_sub15') lat2breakTTL = EEG.event(id2Break+122).latency; end
if strcmp(sub, 'c_sub16') lat2breakTTL = EEG.event(id2Break+169).latency; end
if strcmp(sub, 'c_sub18') lat2breakTTL = EEG.event(id2Break+151).latency; end
if strcmp(sub, 'c_sub19') lat2breakTTL = EEG.event(id2Break+26).latency; end
if strcmp(sub, 'c_sub23') lat2breakTTL = EEG.event(id2Break+71).latency; end
if strcmp(sub, 'c_sub25') lat2breakTTL = EEG.event(id2Break+32).latency; end
if strcmp(sub, 'c_sub29') lat2breakTTL = EEG.event(id2Break+97).latency; end

lat_prev = trlinfo(:, 10:14); 
latencies = reshape (lat_prev', 1, [])'; 
latencies(latencies == 0) = nan; 

sCod = string(trlinfo(:, 1:9));
sCod (ismissing(sCod)) = 'nan';
typOT = ['T' 'V' 'C' 'U' 'R']';
typOT = repmat(typOT, 192, 1);
sCod1 = repelem(sCod, 5, 1); 
%sCod1 = [sCod1, string(typOT)]; %uncomment to include trialinfo
sCodF = join(sCod1, '_');

[d id2Break] = max(abs(diff(latencies)));
lat2breakLog = latencies(id2Break+2)/factor;
if strcmp(sub, 'c_sub09') lat2breakLog = latencies(id2Break+1)/factor; end
newBreak = lat2breakTTL - lat2breakLog

lat2u = 1; 
if strcmp(sub, 'c_sub06') lat2u = 58; end
if strcmp(sub, 'c_sub09') lat2u = 11; end
lat = EEG.event(lat2u).latency; 
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
diff_12 =  latencies(2) / factor - lat;

for i = 1:id2Break
    EEG.event(i) = struct('latency', latencies(i)/factor -diff_12 , 'type', sCodF(i,:), 'urevent', 0);
end
for i = id2Break+1:length(sCodF)
     EEG.event(i) = struct('latency', latencies(i)/factor + newBreak , 'type',sCodF(i,:), 'urevent', 0);
end

%EEG.event = nestedSortStruct(EEG.event, 'latency');
cd(currentPath);



%%

eventChannel = 'POL DC12';if strcmp(sub ,'c_sub09') eventChannel = 'POL DC09'; end
TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);

%% ONE AT THE TIME PARIS

%eeglab
clear, close all
paths = load_paths; 
currentPath = pwd; 

sub = 'p_sub02'

cd(paths.trlinfo)
log_list = dir('*mat'); log_list = {log_list.name};

cd(paths.raw_data)
edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
edf_list = {edf_list.name}';

cd(paths.trlinfo)
load ([sub '_trlinfo.mat']);

cd([paths.raw_data sub '/ieeg/'])
dataset=strcat(paths.raw_data,sub,'/ieeg/');
hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')

data_raw=ft_read_data(dataset);

event=ft_read_event(dataset);

EEG = []; 
EEG.data = data_raw; 
EEG.srate = hdr.Fs;
EEG.chanlocs.label = hdr.label;
EEG.pnts = size(EEG.data,2);
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(event)
    sample =  (event(i).timestamp-double(hdr.FirstTimeStamp))./hdr.TimeStampPerSample + 1;
    EEG.event(i) = struct('latency', sample, 'type', 1, 'urevent', 0);
end

%%
chanids = [1:10];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',[1:10], ...
    'winlength', 100, 'spacing', 10000, 'events', EEG.event);

%%

% check triggers and save in datainfo
% vidA 101, vidB 102, vidC103
ind_vidonset=find([event(:).value]==101|[event(:).value]==102|[event(:).value]==103);
trigger_value=[event(ind_vidonset).value];
trigger_sp=[event(ind_vidonset).sample];

if numel(ind_vidonset)~=size(trlinfo,1)
    error('check why number of trigger do not match number of trials in trlinfo')
end
    
% check diff times between triggers
%trigger_diff=diff(trigger_sp)'.*(1000/EEG.srate); %convert ms
%trlinfo_diff=diff(trlinfo(:,11))./10; % convert ms

%trigger_sp_ds = (trigger_sp).*(1000/EEG.srate); 
trigger_sp_ds = (trigger_sp).*(1000/4096); 
trlinfo_sp_ds = trlinfo(:,11) ./ 10; 

%Downsample data 
%https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
[P,Q] = rat(1000/EEG.srate);
clear data
for chani = 1:size(EEG.data,1)
    chani
    data(chani,:) = resample(double(EEG.data(chani,:)), P, Q);
end

EEG.data = data;
EEG.srate = 1000;
EEG.pnts = size(EEG.data,2);
EEG.xmax = size(EEG.data,2)/EEG.srate;
%EEG.times = 0:2:size(EEG.data,2)*2;

diff_12 = trigger_sp_ds(1) - trlinfo_sp_ds(1); 
% add triggers from trlinfo and event data 
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(trigger_sp)
    EEG.event(i) = struct('latency', trigger_sp_ds(i), 'type', 'trigger', 'urevent', 0);
end
for i = 1:length(trlinfo)
    EEG.event(end+1) = struct('latency',trlinfo_sp_ds(i) + diff_12, 'type', 'trlinfo', 'urevent', 0);
end
EEG.event = nestedSortStruct(EEG.event, 'latency');


%%
chanids = [1:10];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',[1:10], ...
    'winlength', 100, 'spacing', 10000, 'events', EEG.event);

%% LOOP DOWNSAMPLE PARIS

%eeglab
clear, close all
paths = load_paths; 
currentPath = pwd; 


allsubs = {'p_sub01','p_sub02','p_sub03','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14', 'p_sub16','p_sub17'};


for subji = 1:length(allsubs)
        
    clearvars -except paths currentPath allsubs subji

    sub = allsubs{subji};
    
    cd(paths.trlinfo)
    log_list = dir('*mat'); log_list = {log_list.name};
    load ([sub '_trlinfo.mat']);

    cd(paths.raw_data)
    edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
    edf_list = {edf_list.name}';
    
    cd([paths.raw_data sub '/ieeg/'])
    dataset=strcat(paths.raw_data,sub,'/ieeg/');
    hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')
    
    data_raw=ft_read_data(dataset);
    event=ft_read_event(dataset);

    
    EEG = []; 
    EEG.data = data_raw; 
    EEG.srate = hdr.Fs;
    EEG.chanlocs.label = hdr.label;
    EEG.pnts = size(EEG.data,2);
    
    % check triggers and save in datainfo
    % vidA 101, vidB 102, vidC103
    ind_vidonset=find([event(:).value]==101|[event(:).value]==102|[event(:).value]==103);
    trigger_value=[event(ind_vidonset).value];
    trigger_sp=[event(ind_vidonset).sample];

    if strcmp(sub,'p_sub05') ind_vidonset(1:8)=[]; trigger_value(1:8)=[]; trigger_sp(1:8)=[]; end
    
    if numel(ind_vidonset)~=size(trlinfo,1)
        error('check why number of trigger do not match number of trials in trlinfo')
    end
        
    % check diff times between triggers
    %trigger_diff=diff(trigger_sp)'.*(1000/EEG.srate); %convert ms
    %trlinfo_diff=diff(trlinfo(:,11))./10; % convert ms
    
    trigger_sp_ds = (trigger_sp).*(1000/EEG.srate); 
    %trigger_sp_ds = (trigger_sp).*(1000/4096); 
    trlinfo_sp_ds = trlinfo(:,11) ./ 10; 
    
    %Downsample data 
    %https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
    [P,Q] = rat(1000/EEG.srate);
    clear data
    for chani = 1:size(EEG.data,1)
        data(chani,:) = resample(double(EEG.data(chani,:)), P, Q);
    end
    
    EEG.data = data;
    EEG.srate = 1000;
    EEG.pnts = size(EEG.data,2);
    EEG.xmax = size(EEG.data,2)/EEG.srate;
    %EEG.times = 0:2:size(EEG.data,2)*2;
    
    diff_12 = trigger_sp_ds(1) - trlinfo_sp_ds(1); 
    % add triggers from trlinfo and event data 
    EEG.event = struct('latency', [], 'type', '', 'urevent', []);
    for i = 1:length(trigger_sp)
        EEG.event(i) = struct('latency', trigger_sp_ds(i), 'type', 'trigger', 'urevent', 0);
    end
    for i = 1:length(trlinfo)
        EEG.event(end+1) = struct('latency',trlinfo_sp_ds(i) + diff_12, 'type', 'trlinfo', 'urevent', 0);
    end
    EEG.event = nestedSortStruct(EEG.event, 'latency');


    cd ..

    filename = [sub '_dSiEEG.mat']
    save(filename, 'EEG', '-v7.3');
    cd (currentPath)
end






%%
clearvars -except EEG subji paths trlinfo path_log currentPath sub
eventChannel = 'POL DC12';
if strcmp(sub, 'c_sub09') eventChannel = 'POL DC09'; end

TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X > 100000';
EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on');

factor = ( 10000 / EEG.srate); 
x = [EEG.event.latency]';
[d id2Break] = max(abs(diff(x)));

if strcmp(sub, 'c_sub04') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-2);end
if strcmp(sub, 'c_sub06') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-1);end
if strcmp(sub, 'c_sub24') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-4);end
if strcmp(sub, 'c_sub25') [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-8);end

lat2breakTTL = EEG.event(id2Break+1).latency;
if strcmp(sub, 'c_sub09') lat2breakTTL = EEG.event(id2Break+10).latency; end
if strcmp(sub, 'c_sub12') lat2breakTTL = EEG.event(id2Break+63).latency; end
if strcmp(sub, 'c_sub15') lat2breakTTL = EEG.event(id2Break+122).latency; end
if strcmp(sub, 'c_sub16') lat2breakTTL = EEG.event(id2Break+169).latency; end
if strcmp(sub, 'c_sub18') lat2breakTTL = EEG.event(id2Break+151).latency; end
if strcmp(sub, 'c_sub19') lat2breakTTL = EEG.event(id2Break+26).latency; end
if strcmp(sub, 'c_sub23') lat2breakTTL = EEG.event(id2Break+71).latency; end
if strcmp(sub, 'c_sub25') lat2breakTTL = EEG.event(id2Break+32).latency; end
if strcmp(sub, 'c_sub29') lat2breakTTL = EEG.event(id2Break+97).latency; end

lat_prev = trlinfo(:, 10:14); 
latencies = reshape (lat_prev', 1, [])'; 
latencies(latencies == 0) = nan; 

sCod = string(trlinfo(:, 1:9));
sCod (ismissing(sCod)) = 'nan';
typOT = ['T' 'V' 'C' 'U' 'R']';
typOT = repmat(typOT, 192, 1);
sCod1 = repelem(sCod, 5, 1); 
%sCod1 = [sCod1, string(typOT)]; %uncomment to include trialinfo
sCodF = join(sCod1, '_');

[d id2Break] = max(abs(diff(latencies)));
lat2breakLog = latencies(id2Break+2)/factor;
if strcmp(sub, 'c_sub09') lat2breakLog = latencies(id2Break+1)/factor; end
newBreak = lat2breakTTL - lat2breakLog

lat2u = 1; 
if strcmp(sub, 'c_sub06') lat2u = 58; end
if strcmp(sub, 'c_sub09') lat2u = 11; end
lat = EEG.event(lat2u).latency; 
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
diff_12 =  latencies(2) / factor - lat;

for i = 1:id2Break
    EEG.event(i) = struct('latency', latencies(i)/factor -diff_12 , 'type', sCodF(i,:), 'urevent', 0);
end
for i = id2Break+1:length(sCodF)
     EEG.event(i) = struct('latency', latencies(i)/factor + newBreak , 'type',sCodF(i,:), 'urevent', 0);
end

%EEG.event = nestedSortStruct(EEG.event, 'latency');
cd(currentPath);



%%

eventChannel = 'POL DC12';if strcmp(sub ,'c_sub09') eventChannel = 'POL DC09'; end
TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);












%%

