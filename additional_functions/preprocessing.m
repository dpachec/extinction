%% Process EDF files China
%%
%eeglab
clear, close all
paths = load_paths_EXT; 
currentPath = pwd; 
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' }';


for subji = 1 %18:length(allsubs)

    sub = allsubs{subji}; 
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
    
    clearvars -except EEG subji paths trlinfo path_log currentPath sub allsubs EEG1
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
    if strcmp(sub, 'c_sub21') lat2breakTTL = EEG.event(id2Break+122).latency; end
    if strcmp(sub, 'c_sub23') lat2breakTTL = EEG.event(id2Break+71).latency; end
    if strcmp(sub, 'c_sub25') lat2breakTTL = EEG.event(id2Break+32).latency; end    
    if strcmp(sub, 'c_sub26') lat2breakTTL = EEG.event(id2Break+61).latency; end
    if strcmp(sub, 'c_sub29') lat2breakTTL = EEG.event(id2Break+97).latency; end
    if strcmp(sub, 'c_sub30') lat2breakTTL = EEG.event(id2Break+30).latency; end
    
    lat_prev = trlinfo(:, 10:14); 
    latencies = reshape (lat_prev', 1, [])'; 
    latencies(latencies == 0) = nan; 
    
    sCod = string(trlinfo(:, 1:9));
    sCod (ismissing(sCod)) = 'nan';
    typOT = ['T' 'V' 'C' 'U' 'R']';
    typOT = repmat(typOT, 192, 1);
    sCod1 = repelem(sCod, 5, 1); 
    sCod1 = [sCod1, string(typOT)]; %uncomment to include trialinfo
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
    
    
    % % % %resample 
    % % % %  % Downsample data  %https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
    if EEG.srate ~= 1000
        factor = ( 1000 / EEG.srate); 
        [P,Q] = rat(1000/EEG.srate);
        for chani = 1:size(EEG.data, 1)
           chani
           data(chani,:) = resample(double(EEG.data(chani,:)), P, Q); %use the function in fieldtrip, matlab signal processing rescales the data
        end
        events = EEG.event; 
        chans = EEG.chanlocs; 
        EEG = []; 
        EEG.data = data;
        EEG.srate = 1000;
        EEG.pnts = size(EEG.data,2);
        EEG.event = events; 
        x = [events.latency]; 
        x = num2cell(x*factor); 
        [EEG.event.latency] = x{:};
        EEG.chanlocs = chans; 
    end
    
    
    mkdir([paths.ds sub])
    cd ([paths.ds sub])
   
    filename = [sub '_downSampiEEG.mat']
    save(filename, 'EEG', '-v7.3');
    cd (currentPath)
  
    disp('done')
    

    eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
    eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
        'winlength', 50, 'spacing', 10000000, 'events', EEG.event);

   
end


%% Process NCS files PARIS 
%eeglab
clear, close all
paths = load_paths_EXT; 
cd (paths.github)

allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21', 'p_sub22'}'; %subject 8 has differnet format (see below)

% allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
%             'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ... 
%             'p_sub19', 'p_sub20', 'p_sub21'}'; %subject 8 has differnet format (see below)

for subji = 18:18 %1:length(allsubs)

    clearvars -except allsubs subji paths
    sub = allsubs{subji}; 
    cd(paths.trlinfo)
    log_list = dir('*mat'); log_list = {log_list.name};
    
    cd(paths.raw_data)
    edf_list = dir(); edf_list = edf_list(3:end); iDir = [edf_list.isdir]; edf_list = edf_list(iDir);
    edf_list = {edf_list.name}';
    
    cd(paths.trlinfo)
    load ([sub '_trlinfo.mat']);
    
    
    if ~(strcmp(sub, 'p_sub04') | strcmp(sub, 'p_sub15') | strcmp(sub, 'p_sub18'))
        dataset=strcat(paths.raw_data,sub,'/ieeg/');
        hdr=ft_read_header(dataset,'headerformat','neuralynx_ds')
    else
        dataset1 = strcat(paths.raw_data,sub,'/ieeg/sess1/');
        dataset2 = strcat(paths.raw_data,sub,'/ieeg/sess2/');
        hdr=ft_read_header(dataset1,'headerformat','neuralynx_ds')
    end


    %Downsample data  %https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
    %use the function in FS, matlab signal processing rescales the data
    [P,Q] = rat(1000/hdr.Fs);
    nBatch = 12; 
    dataF = zeros(1, (hdr.nSamples * (P/Q)) ); %to test with 1 chan only
    if strcmp(sub, 'p_sub11') dataF = zeros(1, (hdr.nSamples * (P/Q)) -1 ); end %why?
    if strcmp(sub, 'p_sub14') dataF = zeros(1, (hdr.nSamples * (P/Q)) +3 ); end %why?
    if strcmp(sub, 'p_sub04') | strcmp(sub, 'p_sub15') | strcmp(sub, 'p_sub18') dataF = zeros(1, (hdr.nSamples * (P/Q)) *2) ; end 
    eCom = []; 
    for chani = 1:hdr.nChans
        chani
        if ~(strcmp(sub, 'p_sub04') | strcmp(sub, 'p_sub15')| strcmp(sub, 'p_sub18'))
            data_raw_double=ft_read_data(dataset, 'chanindx', chani);
            %data_raw = single(data_raw_double); 
            %clear data_raw_double %save memory
            data_raw = data_raw_double; 
            
        else
            data_raw_double1=ft_read_data(dataset1, 'chanindx', chani);
            data_raw_double2=ft_read_data(dataset2, 'chanindx', chani);
            data_raw = [single(data_raw_double1) single(data_raw_double2)]; 
            clear data_raw_double1 clear data_raw_double2 %save memory %save memory
        end

        hd = length(data_raw) / nBatch; 
        for batchi = 1:nBatch
            if batchi ==1
                data(batchi, :) = resample(data_raw(:,1:hd), P, Q); %use the function in fieldtrip, matlab signal processing rescales the data
            else
                data(batchi, :) = resample(data_raw(:,hd*(batchi-1)+1:hd*batchi), P, Q);
            end
            newC = ['squeeze(data(' num2str(batchi) ', :)) '];
            eCom = [eCom newC]; %build eval command
        end
        
        eval(['dataF = data''; dataF = dataF(:);'])

        EEG = []; 
        EEG.data = dataF'; 
        EEG.srate = 1000; 
        EEG.chanlocs(chani).labels = hdr.label(chani);
        EEG.pnts = size(EEG.data,2);
        EEG.event = struct('latency', [], 'type', '');
    
    if chani == 1 %only store event info in first channel
        if ~(strcmp(sub, 'p_sub04') | strcmp(sub, 'p_sub15') | strcmp(sub, 'p_sub18') )
            event=ft_read_event(dataset);
        end
        if strcmp(sub, 'p_sub18')
            event1=ft_read_event(dataset1); % event.mrk and event.nev must exist in the folder
            event2=ft_read_event(dataset2); % event.mrk and event.nev must exist in the folder
            x = [event2.sample]; x = x + double(hdr.nSamples); 
            x = num2cell(x');
            [event2.sample] = x{:}; 
            event = event1; 
            for i = 1:length(event2)
                event(end+1) = event2(i);
            end
        end
         if strcmp(sub, 'p_sub04') | strcmp(sub, 'p_sub15')
            event=ft_read_event(dataset1); % event.mrk and event.nev must exist in the folder
         end
        

        values = [event.value]; 
        event = event(values~= 0);
        trigger_sp=[event.sample];
        trigger_sp_ds = (trigger_sp).*(1000/hdr.Fs); 
        
    
        EEG.event = struct('latency', [], 'type', []);
        for i = 1:length(trigger_sp)
           % EEG.event(i) = struct('latency', trigger_sp_ds(i), 'type', char(string(event(i).value)));
            EEG.event(i) = struct('latency', trigger_sp_ds(i), 'type', 'trigger');
        end
    
        if strcmp(sub, 'p_sub05') event(1:152) = []; end
        sP1 = find([event.value] == 101); sP1 = sP1(1);
        if strcmp(sub, 'p_sub09') sP1 = find([event.value] == 101);  sP1 = sP1(12); end %experiment restarted
        sP1L = event(sP1).sample*(1000/hdr.Fs);  
        sP2 = find([event.value] == 103); sP2 = sP2(1);
        sP2L = event(sP2).sample*(1000/hdr.Fs);  
        
       
        factor = ( 10000 / EEG.srate); 
        lat_prev = trlinfo(:, 10:14); 
        trlinfo_sp_ds= reshape (lat_prev', 1, [])'; 
        trlinfo_sp_ds(trlinfo_sp_ds == 0) = nan; 
        trlinfo_sp_ds = trlinfo_sp_ds /factor; 
        
        sCod = string(trlinfo(:, 1:9));
        sCod (ismissing(sCod)) = 'nan';
        typOT = ['T' 'V' 'C' 'U' 'R']';
        typOT = repmat(typOT, 192, 1);
        sCod1 = repelem(sCod, 5, 1); 
        sCod1 = [sCod1, string(typOT)]; %uncomment to include trialinfo
        sCodF = join(sCod1, '_');
        sCodF = char(sCodF);
       
        diff_12 = sP1L-trlinfo_sp_ds(2); 
    
        for i = 1:720
            EEG.event(end+1) = struct('latency',  trlinfo_sp_ds(i) + diff_12 , 'type', sCodF(i,:));
            %EEG.event(end+1) = struct('latency', trlinfo_sp_ds(i) + diff_12 , 'type', '2');
        end
        diff_12 = sP2L-trlinfo_sp_ds(722); 
        for i = 721:960
            EEG.event(end+1) = struct('latency',  trlinfo_sp_ds(i) + diff_12 , 'type', sCodF(i,:));
            %EEG.event(end+1) = struct('latency', trlinfo_sp_ds(i) + diff_12 , 'type', '2');
        end
    
        EEG.event = nestedSortStruct(EEG.event, 'latency');

    end



    % % % % % cut data
    t2c = 10000; %10 secs before and after
    trlinfo_sp_ds(isnan(trlinfo_sp_ds)) = []; 
    task_sp = [sP1L - t2c, trlinfo_sp_ds(end) + diff_12  + t2c]; 
    EEG.data = EEG.data(:, task_sp(1):task_sp(2));
    EEG.pnts = size(EEG.data,2);
    x =  [EEG.event.latency]' - ( sP1L - t2c);
    x = num2cell(x');
    EEG.event = struct('latency',x, 'type', {EEG.event.type});
    

    foldN2 = strcat(paths.ds,sub);
    mkdir(foldN2);
    cd(foldN2)
 
     filename = [sub '_' num2str(chani, '%03.f'), '_diEEG.mat'];
     save(filename, 'EEG', '-v7.3');
     cd (paths.github);

    end

% % 
% % % % % % % only makes sense for the first channel
% % chanids = [1];
% % eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',chanids, ...
% %     'winlength', 50, 'spacing', 10000, 'events', EEG.event);




end



%% just 2 check, plot with original SR (in 4096) and with no cuts

clear, close all
paths = load_paths_EXT; 
cd (paths.github)

allsubs = {'p_sub01','p_sub02','p_sub03','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14', 'p_sub16','p_sub17'};


sub = allsubs{1}
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

for chani = 1:1
    data_raw_double=ft_read_data(dataset, 'chanindx', chani);
    data_raw = single(data_raw_double); 
    clear data_row_double %save memory
end


EEG = []; 
EEG.data = data_raw; 
EEG.srate = hdr.Fs; 
EEG.chanlocs.labels = hdr.label;
EEG.pnts = size(EEG.data,2);
EEG.event = struct('latency', [], 'type', '');

event=ft_read_event(dataset);
values = [event.value]; 
event = event(values~= 0);
trigger_sp=[event.sample];

EEG.event = struct('latency', [], 'type', []);
for i = 1:length(trigger_sp)
   % EEG.event(i) = struct('latency', trigger_sp(i), 'type', char(string(event(i).value)));
    EEG.event(i) = struct('latency', trigger_sp(i), 'type', 'trigger');
end




%% IMPORTANT > group the multiple electrode file in one channel (Paris data)

clear, close all
paths = load_paths_EXT; 
cd (paths.github)

allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ...
            'p_sub19', 'p_sub20', 'p_sub21', 'p_sub22'}'; %subject 8 has differnet format (see below)



for subji = 18:20 %:length(allsubs)

    clearvars -except allsubs subji paths
    sub = allsubs{subji}; 
    
    cd([paths.ds sub ])
    diEEG_list = dir(); diEEG_list = diEEG_list(3:end); diEEG_list = {diEEG_list.name}';
    
    for chani = 1:length(diEEG_list)
        load(diEEG_list{chani})
        if chani ==1
            events = EEG.event
        end
        data(chani,:) = EEG.data; 
        chans2rec(chani) = EEG.chanlocs(chani).labels;

    end
    
    EEG.data = data; 
    [EEG.chanlocs.labels] = chans2rec{:}; 
    EEG.event = events;

    filename = [sub '_downSampiEEG.mat']
    save(filename, 'EEG', '-v7.3');


end




%% COMBINE ELECTRODE AND IEEG DATA, FILTER AND REMOVE ARTIFACTS (CRITICAL STEP)

clear, close all
paths = load_paths_EXT; 
currentPath = paths.github; 
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ...
            'p_sub19', 'p_sub20', 'p_sub21', 'p_sub22'}'; %subject 8 has differnet format (see below)}';


for subji = 50:50 % 1:length(allsubs)

    clearvars -except allsubs subji paths
    sub = allsubs{subji}; 
    cd([paths.ds sub ])
    diEEG_list = dir('*downSampiEEG.mat');
    diEEG_list = {diEEG_list.name}';
    load(diEEG_list{1})
    %load electrode files
    cd([paths.raw_data sub '\elec'])


    if sub(1) == 'c' %guangzhou data
        chanL= dir('*.mat');chanL= {chanL.name}';
        load(chanL{1})
        elec_info = table2cell(elec_info);
        for chani = 1:size(elec_info, 1)
            
            chanLabel = elec_info(chani, 1);
            chanEEG = {EEG.chanlocs.labels}'; 
            chanEEG = erase(chanEEG, 'POL');
            chanEEG = erase(chanEEG, ' ');
            chanMatch = strmatch(chanLabel, chanEEG, 'exact');
            if length(chanMatch) ==1
                if ~isempty(chanMatch)
                    EEG.chanlocs(chanMatch).fsLabel = elec_info{chani,2};
                    EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,9});
                    if subji == 17
                        EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,10});
                    end
                end
            end
        end
    end
    
    
    if sub(1) == 'p' %paris Data
        chanL= dir('*.csv');chanL= {chanL.name}';
        elec_info = readtable(chanL{1}, 'Delimiter', ';')
        elec_info1 = elec_info; 
        elec_info1(:, 1) = fillmissing(elec_info(:, 1),'previous')
        elec_info1 = table2cell(elec_info1);
        for chani = 1:size(elec_info1, 1)    
            if ~isnan(elec_info1{chani, 2})
                chanLabel = strcat(elec_info1(chani, 1), '_', string(elec_info1(chani, 2)));
                % % % correct should be labels and not label, this field name was incorrect in the previous preprocessing stage (important
                % for epoching data later)
                if isfield(EEG.chanlocs, 'label')
                    chanEEG = {EEG.chanlocs.label}'; 
                    [EEG.chanlocs.labels] = EEG.chanlocs.label;
                    EEG.chanlocs = rmfield(EEG.chanlocs,'label');
                else 
                    chanEEG = {EEG.chanlocs.labels}';
                end
        
                chanMatch = strmatch(chanLabel, chanEEG, 'exact'); 
                if ~isempty(chanMatch)
                    EEG.chanlocs(chanMatch).fsLabel = elec_info1{chani,36};
                    EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,4:6});
                end
            end
        end
    end

    % % % % re-reference 
    EEG = re_reference_EXT(EEG, 'bipo', sub);

    EEG = manual_channel_removal(EEG, sub); 
    
    % % %  Filter
    EEG = add_EEGLAB_fields(EEG);
    EEG = pop_eegfiltnew (EEG, .1, 300);
    EEG = pop_eegfiltnew (EEG, 49, 50, [],  1); %notch filter
    EEG = pop_eegfiltnew (EEG, 99, 101, [],  1); %notch filter
    EEG = pop_eegfiltnew (EEG, 149, 150, [],  1); %notch filter

    % % % % remove artifacts
    EEG = artifact_detection_EXT(EEG, 6, 4, 1000, 500);

   
    % % % % save final versions
    mkdir(paths.iEEG)
    cd([paths.iEEG])
    filename = [sub '_iEEG.mat']
    %save(filename, 'EEG', '-v7.3');
    save(filename, 'EEG');


end


%% EXPORT ELECTROEDS JUST TO TEST

clear, close all
paths = load_paths_EXT; 
currentPath = paths.github; 
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


for subji = 1:length(allsubs)

    clearvars -except allsubs subji paths ALLEEG
    sub = allsubs{subji}; 
    cd([paths.ds sub ])
    diEEG_list = dir('*downSampiEEG.mat');diEEG_list = {diEEG_list.name}';
    load(diEEG_list{1})
    %load electrode files
    cd([paths.raw_data sub '\elec'])


    if sub(1) == 'c' %guangzhou data
        chanL= dir('*.mat');chanL= {chanL.name}';
        load(chanL{1})
        elec_info = table2cell(elec_info);
        for chani = 1:size(elec_info, 1)
            
            chanLabel = elec_info(chani, 1);
            chanEEG = {EEG.chanlocs.labels}'; 
            chanEEG = erase(chanEEG, 'POL');
            chanEEG = erase(chanEEG, ' ');
            chanMatch = strmatch(chanLabel, chanEEG, 'exact');
            if length(chanMatch) ==1
                if ~isempty(chanMatch)
                    EEG.chanlocs(chanMatch).fsLabel = elec_info{chani,2};
                    EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,9});
                    if subji == 17
                        EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,10});
                    end
                end
            end
        end
    end
    
    
    if sub(1) == 'p' %paris Data
        chanL= dir('*.csv');chanL= {chanL.name}';
        elec_info = readtable(chanL{1}, 'Delimiter', ';')
        elec_info1 = elec_info; 
        elec_info1(:, 1) = fillmissing(elec_info(:, 1),'previous')
        elec_info1 = table2cell(elec_info1)
        for chani = 1:size(elec_info1, 1)    
            if ~isnan(elec_info1{chani, 2})
                chanLabel = strcat(elec_info1(chani, 1), '_', string(elec_info1(chani, 2)));
                % % % correct should be labels and not label, this field name was incorrect in the previous preprocessing stage (important
                % for epoching data later)
                if isfield(EEG.chanlocs, 'label')
                    chanEEG = {EEG.chanlocs.label}'; 
                    [EEG.chanlocs.labels] = EEG.chanlocs.label;
                    EEG.chanlocs = rmfield(EEG.chanlocs,'label');
                else 
                    chanEEG = {EEG.chanlocs.labels}';
                end
        
                chanMatch = strmatch(chanLabel, chanEEG, 'exact'); 
                if ~isempty(chanMatch)
                    EEG.chanlocs(chanMatch).fsLabel = elec_info1{chani,36};
                    EEG.chanlocs(chanMatch).mniCoord = double(elec_info{chani,4:6});
                end
            end
        end
    end

    % % % % re-reference 
    EEG = re_reference_EXT(EEG, 'bipo', sub);

    EEG = manual_channel_removal(EEG, sub);
    
   
   
    EEG = rem_EEGLAB_fields(EEG);

    ALLEEG{subji, :} = EEG; 


end


















%% save all electrodes into a csv file

clear, close all
paths = load_paths; 

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


allSChans = []; 
for subji = 1:length(allsubs)
    
    sub = allsubs{subji}; 
    cd([ paths.fiEEG])
    load ([sub '_iEEG.mat']);
    
    chansLab = struct2cell(EEG.chanlocs)';
    chansLab(:, 4) = {sub}; %{num2str(subji, '%02.f')};

    %allSChans = [allSChans; chansLab]; 

    allSChans{subji} = chansLab; 
end

cd (paths.elec)
save('allSChans', 'allSChans')








%% subject 8 in paris has a different format MICROMED % discuss data inconsistencies with Katia

clear, close all
paths = load_paths; 

sub = 'p_sub08'

hdr = ft_read_header('D:\extinction\raw_data\p_sub08\ieeg\EEG_444702.TRC');
data = ft_read_data('D:\extinction\raw_data\p_sub08\ieeg\EEG_444702.TRC');
EEG = []; 
EEG.data = data;
EEG.srate = hdr.Fs;
EEG.pnts = size(EEG.data,2);
EEG.chanlocs = struct('labels', hdr.label)
EEG = add_EEGLAB_fields(EEG);


 
cd(paths.trlinfo)
load ([sub '_trlinfo.mat']);


%%

clearvars -except EEG subji paths trlinfo path_log currentPath sub allsubs EEG1
eventChannel = 'sti1+';
TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X > 180 & X < 400';
EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off','delevent', 'on');

chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 1000, 'events', EEG.event);


%%


EEG.event = struct('latency', [], 'type', '', 'urevent', []);

factor = ( 10000 / EEG.srate); 
x = [EEG.event.latency]';
[d id2Break] = max(abs(diff(x)));

lat_prev = trlinfo(:, 10:14); 
latencies = reshape (lat_prev', 1, [])'; 
latencies(latencies == 0) = nan; 

sCod = string(trlinfo(:, 1:9));
sCod (ismissing(sCod)) = 'nan';
typOT = ['T' 'V' 'C' 'U' 'R']';
typOT = repmat(typOT, 192, 1);
sCod1 = repelem(sCod, 5, 1); 
sCod1 = [sCod1, string(typOT)]; %uncomment to include trialinfo
sCodF = join(sCod1, '_');

[d id2Break] = max(abs(diff(latencies)));
lat2breakLog = latencies(id2Break+2)/factor;

diff_12 = 5150100;
%diff_13 = 6721200; 
for i = 1:380
    EEG.event(i) = struct('latency', latencies(i)/factor  + diff_12 , 'type', sCodF(i,:), 'urevent', 0);
end
for i = 381:395
    EEG.event(i) = struct('latency', latencies(i)/factor  + diff_12-1700 , 'type', sCodF(i,:), 'urevent', 0);
end
for i = 396:488
    EEG.event(i) = struct('latency', latencies(i)/factor  + diff_12 -3500, 'type', sCodF(i,:), 'urevent', 0);
end
for i = 489:720
    EEG.event(i) = struct('latency', latencies(i)/factor  + diff_12 -7000, 'type', sCodF(i,:), 'urevent', 0);
end

% for i = 721:length(sCodF)
%      EEG.event(i) = struct('latency', latencies(i)/factor + diff_13, 'type',sCodF(i,:), 'urevent', 0);
% end


eventChannel = 'sti1+'; TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 1000, 'events', EEG.event);




%%
% % % %  resample 
%Downsample data  %https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
if EEG.srate ~= 1000
    factor = ( 1000 / EEG.srate); 
    [P,Q] = rat(1000/EEG.srate);
    for chani = 1:size(EEG.data, 1)
       chani
       data(chani,:) = resample(double(EEG.data(chani,:)), P, Q); %use the function in fieldtrip, matlab signal processing rescales the data
    end
    events = EEG.event; 
    chans = EEG.chanlocs; 
    EEG = []; 
    EEG.data = data;
    EEG.srate = 1000;
    EEG.pnts = size(EEG.data,2);
    EEG.event = events; 
    x = [events.latency]; 
    x = num2cell(x*factor); 
    [EEG.event.latency] = x{:};
    EEG.chanlocs = chans; 
end


mkdir([paths.ds sub])
cd ([paths.ds sub])

filename = [sub '_downSampiEEG.mat']
save(filename, 'EEG', '-v7.3');
cd (currentPath)

disp('done')


eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);

%%

factor = ( 1000 / hdr.Fs); 
[P,Q] = rat(factor);
for chani = 1:size(data, 1)
   chani
   dataRS(chani,:) = resample(double(data(chani,:)), P, Q); %use the function in fieldtrip, matlab signal processing rescales the data
end
events = EEG.event; 
chans = EEG.chanlocs; 







%% plot example

eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];

eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);






















































