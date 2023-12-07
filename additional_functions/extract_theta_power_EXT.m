 
function [EEG] = extract_theta_power_EXT (EEG)
    
    
    for chani = 1:size(EEG.data, 1)
        for triali = 1:size(EEG.data, 3)

            data = EEG.data(chani, :, triali);

            if isempty(find(isnan(data)))
                [ phaseData, powerData, amplitudeData ]= hilbertTransformation_hui( data, 1000, 3,8);
                EEG.power(triali, chani, 1, :) = powerData ;
            else
                EEG.power(triali, chani, 1, :) = nan(1, length(data));
            end

        end
    end


end

