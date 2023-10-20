function [EEG] = normalize_EXT(EEG)

if isfield(EEG, 'power')
  if ndims(EEG.power) == 3 
      tmp(:, 1, :, :) = EEG.power; 
      EEG.power = tmp; 
  end
  % % % % % % across trials
  mT = mean(EEG.power,1, 'omitnan');
  stdT = std(EEG.power,[], 1, 'omitnan');
  EEG.power = bsxfun(@rdivide, EEG.power- mT, stdT);  
end

if isfield(EEG, 'data')
  % % % % % % across trials
  mT = mean(EEG.data,3, 'omitnan');
  stdT = std(EEG.data,[], 3, 'omitnan');
  EEG.data = bsxfun(@rdivide, EEG.data- mT, stdT);  
end
    

 
end
























                