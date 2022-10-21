function[EEG] = rem_EEGLAB_fields(EEG)

EEG = rmfield(EEG,'setname');
EEG = rmfield(EEG,'filename');
EEG = rmfield(EEG,'filepath');
EEG = rmfield(EEG,'subject');
EEG = rmfield(EEG,'group');
EEG = rmfield(EEG,'session');
EEG = rmfield(EEG,'comments');
EEG = rmfield(EEG,'nbchan');
EEG = rmfield(EEG,'trials');
EEG = rmfield(EEG,'pnts');
EEG = rmfield(EEG,'xmin');
EEG = rmfield(EEG,'xmax');
EEG = rmfield(EEG,'times');
EEG = rmfield(EEG,'icaact');
EEG = rmfield(EEG,'icawinv');
EEG = rmfield(EEG,'icasphere');
EEG = rmfield(EEG,'icaweights');
EEG = rmfield(EEG,'icachansind');
EEG = rmfield(EEG,'urchanlocs');
EEG = rmfield(EEG,'chaninfo');
EEG = rmfield(EEG,'ref');
EEG = rmfield(EEG,'urevent');
EEG = rmfield(EEG,'eventdescription');
EEG = rmfield(EEG,'epoch');
EEG = rmfield(EEG,'epochdescription');
EEG = rmfield(EEG,'condition');
EEG = rmfield(EEG,'reject');
EEG = rmfield(EEG,'stats');
EEG = rmfield(EEG,'specdata');
EEG = rmfield(EEG,'specicaact');
EEG = rmfield(EEG,'splinefile');
EEG = rmfield(EEG,'icasplinefile');
EEG = rmfield(EEG,'dipfit');
EEG = rmfield(EEG,'history');
EEG = rmfield(EEG,'saved');
EEG = rmfield(EEG,'etc');
EEG = rmfield(EEG,'run');

% can also be done by passing all in an array
%EEG.chanlocs = rmfield(EEG.chanlocs, {'ref', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', 'type', 'urchan'});


end