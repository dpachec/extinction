%% % % load edf and trial info file

eeglab
clear, close all
load D:/extinction/data/preproc/trialinfo/c_sub01_trlinfo.mat
%load D:/extinction/data/preproc/trialinfo/c_sub02_trlinfo.mat
%load D:/extinction/data/preproc/trialinfo/c_sub03_trlinfo.mat

EEG = pop_biosig('D:\extinction\raw_data\c_sub01\ieeg\DBX-Learning test.edf', 'importevent', 'off'); 
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
%plot(diff(x))
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


%%



%%

EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(x)
    EEG.event(length(EEG.event) + 1) = struct('latency', x - diff_13, 'type', '1', 'urevent', 0);
end
EEG.event(1) = [];


eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);










%%
figure
plot (latencies)
[d id] = max(abs(diff(latencies)));


%% 
x = [EEG.event.latency]';
figure
plot (diff(x))
[d id] = max(abs(diff(x)));


x(id+1:end) = x(id+1:end) + d + 83000; 


%%

EEG.event = []; 
EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(x)
    EEG.event(length(EEG.event) + 1) = struct('latency', x(i)   , 'type', '1', 'urevent', 0);
end
EEG.event(1) = [];


eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);







%% ALL IN LOOP

eeglab
clear, close all
path_log = 'D:/extinction/data/preproc/trialinfo/'
path_edf = 'D:\extinction\raw_data\'
currentPath = pwd; 

for subji = 1:1

    cd(path_log)
    log_list = dir('*mat'); log_list = {log_list.name};
    load (log_list{subji});
    lat_prev = trlinfo(:, 10:14); 
    latencies = reshape (lat_prev', 1, [])'; 
    latencies = latencies(~isnan(latencies));
    
    cd([path_edf 'c_sub0' num2str(subji) '\ieeg\'])
    edf_list = dir('*edf'); edf_list = {edf_list.name};
    EEG = pop_biosig(edf_list{1}, 'importevent', 'off'); 

    eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
    strEvent = 'X > 100000';
    EEG.event = [];
    EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
        'delevent', 'on');

    
    factor = ( 10000 / EEG.srate); 
    x = [EEG.event.latency]';
    [d id2Break] = max(abs(diff(x)));
    lat2breakTTL = EEG.event(id2Break+1).latency; 
    
    x = latencies;
    [d id2Break] = max(abs(diff(x)));
    %lat2breakLog = 1444000; %x(id) ; 
    lat2breakLog = x(id2Break+2)/factor;
    newBreak = lat2breakTTL - lat2breakLog

    lat = EEG.event(1).latency; 

    EEG.event = struct('latency', [], 'type', '', 'urevent', []);
    diff_12 =  latencies(2) / factor- lat;
    for i = 1:id2Break
        EEG.event(i) = struct('latency', latencies(i)/factor -diff_12 , 'type', '1', 'urevent', 0);
    end
    for i = id2Break+1:length(latencies)
        EEG.event(i) = struct('latency', latencies(i)/factor + newBreak , 'type', '2', 'urevent', 0);
    end

    cd(currentPath);

end

disp ('done');
%%

eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);











