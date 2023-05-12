function[EEG] = manual_channel_removal(EEG, sub)


    if strcmp(sub, 'c_sub06')
        chans2rem = {'POL A''2', 'POL A''3'}; 
        id2rem = contains({EEG.chanlocs.labels}', chans2rem)
        EEG.chanlocs(id2rem) = []; 
        EEG.data(id2rem, :) = []; 
    end

    if strcmp(sub, 'c_sub02')
        chans2rem = {'POL H''8'}; 
        id2rem = contains({EEG.chanlocs.labels}', chans2rem)
        EEG.chanlocs(id2rem) = []; 
        EEG.data(id2rem, :) = []; 
    end
    




end