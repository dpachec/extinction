function [EEG]  = artifact_detection_EXT(EEG, std_thres, std_thres2, paddingValue, minSegLength)

    for chani = 1:size(EEG.data, 1)
    
            data = EEG.data(chani, :);
            zScoreAmp = bsxfun(@rdivide, data - mean(data), std(data));
            grad = diff(data); grad(end+1) = nan; 
            zScoreGrad = bsxfun(@rdivide, grad - mean(grad, 'omitnan'), std(grad, 'omitnan'));
            hpf_data = eegfilt(data,1000, 250, nan);
            zScoreHPFD = bsxfun(@rdivide, hpf_data - mean(hpf_data), std(hpf_data));
            
            markers = zeros(1, length(data));
            
            markers(zScoreAmp > std_thres | zScoreGrad > std_thres | zScoreHPFD > std_thres | ... 
                (zScoreAmp > std_thres2 & (zScoreGrad > std_thres2 | zScoreHPFD > std_thres2)) ) = 1; %check Staresina 2018 Nat Neuro
            
            newTrace = padding_EXT(markers, paddingValue); 
            
            markers = remove_small_segments_EXT(newTrace, minSegLength); 
    
            EEG.data(chani, markers == 1) = nan;
            %EEG.markers_artifacts(chani, :) = markers; 
    
            
            
    end

end











