function [data] = padding_EXT(data, paddingValue)

    clustinfo = bwconncomp(data, 8); %see help bwconncomp for 8
    
    if ~isempty(clustinfo.PixelIdxList) %check origin and end
        for i=1:length(clustinfo.PixelIdxList)
            clust = clustinfo.PixelIdxList{i};
            if clust(1) - paddingValue > 1 & clust(end)+ paddingValue < length(data) 
                data(clust(1)-paddingValue:clust(1)) = 1; 
                data(clust(end):clust(end)+paddingValue) = 1; 
            end
        end
    end

end