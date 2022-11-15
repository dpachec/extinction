 
function [EEG] = extract_power_NAV (EEG, timeRes)
    
    data_ft = mat2ft(EEG.data, EEG.srate); %1000 sr
    

    %fieldtrip analysis
    cfg              = [];
    cfg.method       = 'wavelet'; %%'mtmconvol'; % 'wavelet'; %
    cfg.width        = linspace(3, 6, 29);
    cfg.output       = 'pow';
    cfg.foi          = [1:1:29];   
    cfg.pad          = 'nextpow2';
    if strcmp (timeRes, 'all')
        cfg.toi          = 'all'; % takes as reference the number of time windows defined above
    else
        cfg.toi          = 0.005:timeRes:(length(EEG.data)/1000) - 0.005;
    end
    cfg.keeptrials   = 'yes'; % keep individual trials or average
    cfg.showcallinfo = 'no';% no log console
    cfg.feedback     = 'no'; 
    tf_data_L          = ft_freqanalysis(cfg, data_ft);
    dataL = tf_data_L.powspctrm;


 
    cfg              = [];
    cfg.method       = 'wavelet'; %'mtmconvol'; % 'wavelet'; %
    cfg.width        = linspace(6, 12, 25);
    cfg.output       = 'pow';
    cfg.foi          = [30:5:150];   % analysis 2 to X Hz in steps of 2 Hz 
    cfg.pad          = 'nextpow2';
    if strcmp (timeRes, 'all')
        cfg.toi          = 'all'; % takes as reference the number of time windows defined above
    else
        cfg.toi          = 0.005:timeRes:(length(EEG.data)/1000) - 0.005;
    end
    cfg.showcallinfo = 'no';% no log console
    cfg.keeptrials   = 'yes'; % keep individual trials, if not, it makes an average 
    cfg.feedback     = 'no'; 
    tf_data_H          = ft_freqanalysis(cfg, data_ft);
    dataH = tf_data_H.powspctrm;

    dataLH = cat (3, dataL, dataH);


    EEG.power = squeeze(dataLH);



end

