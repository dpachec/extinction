function [EEG] = extract_phases_EXT(EEG, f)

data = EEG.data;
[phaseData, powerData, amplitudeData] = hilbertTransformation_hui(data, 1000, 3, 8)




% 
% oneListTraces = permute (EEG.data, [2 1 3]); 
% data = mat2ft(oneListTraces, 1000); %1000 SR
% 
% cfg = [];
% cfg.output     = 'all';
% cfg.bpfilter   = 'yes';
% cfg.bpfreq     = f;
% cfg.hilbert    = 'angle';
% cfg.keeptrials = 'yes';
% data_hilbert = ft_preprocessing(cfg,data);
% 
% 
% 
% 
% 
% EEG.phases = cat(3,data_hilbert.trial{:});

