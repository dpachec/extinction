function [EEG ]  = artifact_detection(EEG, nIQR, paddingValue, minSegLength)

countCh = 1;
sr = 1000;

EEG.nbchan = size(EEG.data,1);
markers = zeros (size(EEG.data,1), size(EEG.data,2));

for chani = 1:size(EEG.data, 1)

    dataTmp    = EEG.data(chani,:);
    %dataTmpG   = diff(dataTmp);dataTmpG(end+1) = 0;
    T  = median(dataTmp, 'omitnan') + iqr(dataTmp) * nIQR;
    %TG = median(dataTmpG, 'omitnan') + iqr(dataTmpG) * nIQR;
    
    %EEG.thresholds(chani) = T; %keep the threshold value for each channel
    dataTmp(abs(dataTmp)  > T) = NaN;
    %dataTmp(abs(dataTmpG) > TG) = NaN;

    markersNoise = zeros (1, length(dataTmp)); 
    markersNoise(isnan(dataTmp)) = 1;
    
    newTrace = padding_NAV(markersNoise, paddingValue); 
    
    newTrace = remove_small_segments_NAV(newTrace, minSegLength); 
      
    markers(chani, :) = newTrace;

    EEG.data(chani, newTrace==1) = nan; 
       
end

EEG.marker_artifacts = markers;


end













