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
path_log = '/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/preproc/trialinfo/'
path_edf = '/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/raw_data/'
currentPath = pwd; 

sub = 'c_sub23'

cd(path_log)
log_list = dir('*mat'); log_list = {log_list.name};

cd(path_edf)
edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
edf_list = {edf_list.name}';

cd(path_log)
load ([sub '_trlinfo.mat']);

cd([path_edf sub '/ieeg/'])
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

lat2breakTTL = EEG.event(id2Break+1).latency;
if strcmp(sub, 'c_sub09') lat2breakTTL = EEG.event(id2Break+10).latency; end
if strcmp(sub, 'c_sub12') lat2breakTTL = EEG.event(id2Break+63).latency; end
if strcmp(sub, 'c_sub15') lat2breakTTL = EEG.event(id2Break+122).latency; end
if strcmp(sub, 'c_sub16') lat2breakTTL = EEG.event(id2Break+169).latency; end
if strcmp(sub, 'c_sub18') lat2breakTTL = EEG.event(id2Break+151).latency; end
if strcmp(sub, 'c_sub19') lat2breakTTL = EEG.event(id2Break+26).latency; end

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



%% ALL IN LOOP

%eeglab
clear, close all
path_log = '/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/preproc/trialinfo/'
path_edf = '/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/raw_data/'
currentPath = pwd; 

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub29' };

cd(path_log)
log_list = dir('*mat'); log_list = {log_list.name};

cd(path_edf)
edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
edf_list = {edf_list.name}';
       
for subji = 1:length(allsubs)

    sub = allsubs{subji};
    cd(path_log)
    load ([sub '_trlinfo.mat']);
    
    cd([path_edf sub '/ieeg/'])
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
    
    lat2breakTTL = EEG.event(id2Break+1).latency;
    if strcmp(sub, 'c_sub09') lat2breakTTL = EEG.event(id2Break+10).latency; end
    if strcmp(sub, 'c_sub12') lat2breakTTL = EEG.event(id2Break+63).latency; end
    if strcmp(sub, 'c_sub15') lat2breakTTL = EEG.event(id2Break+122).latency; end
    if strcmp(sub, 'c_sub16') lat2breakTTL = EEG.event(id2Break+169).latency; end
    if strcmp(sub, 'c_sub18') lat2breakTTL = EEG.event(id2Break+151).latency; end
    if strcmp(sub, 'c_sub19') lat2breakTTL = EEG.event(id2Break+26).latency; end
    
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

end

disp ('done');


%%

eventChannel = 'POL DC12';if strcmp(sub ,'c_sub09') eventChannel = 'POL DC09'; end
TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);





%%

p2u = 1

clearvars -except EEG subji paths trlinfo path_log p2u
    eventChannel = 'POL DC12';
    if subji == 9 eventChannel = 'POL DC09'; end
    TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
    strEvent = 'X > 100000';
    EEG.event = [];
    EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
        'delevent', 'on');

    factor = ( 10000 / EEG.srate); 
    x = [EEG.event.latency]';
    [d id2Break] = max(abs(diff(x)));
    
    if subji == 4 [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-2);end
    if subji == 6 [ii,ii] = sort(abs(diff(x))); id2Break = ii(end-1);end
    
    lat2breakTTL = EEG.event(id2Break+1).latency; 
    if subji == 9 lat2breakTTL = EEG.event(id2Break+10).latency; end
    if subji == 12 lat2breakTTL = EEG.event(id2Break+63).latency; end
    if subji == 15 lat2breakTTL = EEG.event(id2Break+122).latency; end
    if subji == 16 lat2breakTTL = EEG.event(id2Break+169).latency; end
    if subji == 18 lat2breakTTL = EEG.event(id2Break+151).latency; end
    if subji == 19 lat2breakTTL = EEG.event(id2Break+26).latency; end
    
    
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
    if subji == 9 lat2breakLog = latencies(id2Break+1)/factor; end
    newBreak = lat2breakTTL - lat2breakLog

    lat = EEG.event(p2u).latency; 
    EEG.event = struct('latency', [], 'type', '', 'urevent', []);
    diff_12 =  latencies(2) / factor - lat;

    for i = 1:id2Break
        EEG.event(i) = struct('latency', latencies(i)/factor -diff_12 , 'type', sCodF(i,:), 'urevent', 0);
    end
    for i = id2Break+1:length(sCodF)
         EEG.event(i) = struct('latency', latencies(i)/factor + newBreak , 'type',sCodF(i,:), 'urevent', 0);
    end
    




eventChannel = 'POL DC12';
if subji == 9 eventChannel = 'POL DC09'; end
TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);





%%

