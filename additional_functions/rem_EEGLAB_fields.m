function[EEG] = rem_EEGLAB_fields(EEG)

%disp('hola')
if isfield (EEG, 'chanlocs') & isfield (EEG.chanlocs, 'ref')
    EEG.chanlocs = rmfield(EEG.chanlocs,{'ref', 'theta', 'radius', 'X', 'Y', 'Z', 'sph_theta', 'sph_phi', 'sph_radius', 'type', 'urchan'});
end

if isfield(EEG, 'setname')
    EEG = rmfield(EEG,'setname');
end
if isfield(EEG, 'filename')
    EEG = rmfield(EEG,'filename');
end
if isfield(EEG, 'filepath')
    EEG = rmfield(EEG,'filepath');
end
if isfield(EEG, 'subject')
    EEG = rmfield(EEG,'subject');
end
if isfield(EEG, 'group')
    EEG = rmfield(EEG,'group');
end
if isfield(EEG, 'session')
EEG = rmfield(EEG,'session');
end

if isfield(EEG, 'comments')
EEG = rmfield(EEG,'comments');
end

if isfield(EEG, 'nbchan')
EEG = rmfield(EEG,'nbchan');
end

if isfield(EEG, 'trials')
EEG = rmfield(EEG,'trials');
end

if isfield(EEG, 'pnts')
EEG = rmfield(EEG,'pnts');
end

if isfield(EEG, 'xmin')
EEG = rmfield(EEG,'xmin');
end

if isfield(EEG, 'xmax')
EEG = rmfield(EEG,'xmax');
end

if isfield(EEG, 'times')
EEG = rmfield(EEG,'times');
end

if isfield(EEG, 'icaact')
EEG = rmfield(EEG,'icaact');
end

if isfield(EEG, 'icawinv')
EEG = rmfield(EEG,'icawinv');
end

if isfield(EEG, 'icasphere')
EEG = rmfield(EEG,'icasphere');
end

if isfield(EEG, 'icaweights')
EEG = rmfield(EEG,'icaweights');
end

if isfield(EEG, 'icachansind')
EEG = rmfield(EEG,'icachansind');
end

if isfield(EEG, 'urchanlocs')
EEG = rmfield(EEG,'urchanlocs');
end

if isfield(EEG, 'chaninfo')
EEG = rmfield(EEG,'chaninfo');
end

if isfield(EEG, 'ref')
EEG = rmfield(EEG,'ref');
end

if isfield(EEG, 'urevent')
EEG = rmfield(EEG,'urevent');
end

if isfield(EEG, 'eventdescription')
EEG = rmfield(EEG,'eventdescription');
end

if isfield(EEG, 'epoch')
EEG = rmfield(EEG,'epoch');
end

if isfield(EEG, 'epochdescription')
EEG = rmfield(EEG,'epochdescription');
end

if isfield(EEG, 'condition')
EEG = rmfield(EEG,'condition');
end

if isfield(EEG, 'reject')
EEG = rmfield(EEG,'reject');
end

if isfield(EEG, 'stats')
EEG = rmfield(EEG,'stats');
end

if isfield(EEG, 'specdata')
EEG = rmfield(EEG,'specdata');
end

if isfield(EEG, 'specicaact')
EEG = rmfield(EEG,'specicaact');
end

if isfield(EEG, 'splinefile')
EEG = rmfield(EEG,'splinefile');
end

if isfield(EEG, 'icasplinefile')
EEG = rmfield(EEG,'icasplinefile');
end

if isfield(EEG, 'dipfit')
EEG = rmfield(EEG,'dipfit');
end

if isfield(EEG, 'history')
EEG = rmfield(EEG,'history');
end

if isfield(EEG, 'saved')
EEG = rmfield(EEG,'saved');
end

if isfield(EEG, 'etc')
EEG = rmfield(EEG,'etc');
end

if isfield(EEG, 'run')
EEG = rmfield(EEG,'run');
end

% can also be done by passing all in an array
%EEG.chanlocs = rmfield(EEG.chanlocs, {'ref', 'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius', 'type', 'urchan'});


end