function [EEG] = normalize_EXT(EEG)

  if ndims(EEG.power) == 3 
      tmp(:, 1, :, :) = EEG.power; 
      EEG.power = tmp; 
  end


  % % % % % % across trials
  mT = mean(EEG.power,1, 'omitnan');
  stdT = std(EEG.power,[], 1, 'omitnan');
  EEG.power = bsxfun(@rdivide, EEG.power- mT, stdT);  


    



%     for triali = 1:size(EEG.power, 1)
%         for chani = 1:size(EEG.power, 2)
% 
%              data = EEG.power(triali, chani, :, :); 
%              mT = mean(data(:, :, :, 151:250),4, 'omitnan');
%              stdT = std(data(:, :, :, 151:250),[], 4, 'omitnan');
%              dataNorm = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);  
% 
% 
%             EEG.power(triali, chani, :, :) = dataNorm; 
% 
%         end
%     end


 
end
























                