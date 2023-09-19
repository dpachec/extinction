function RDM_TS=remDiagRDM_TS(RDM_TS)


    for timei = 1:size(RDM_TS, 3)
        RDM_TS_t = RDM_TS(:, :, timei); 
        RDM_TS_t(eye(size(RDM_TS_t))==1) = nan;
        RDM_TS(:, :, timei) = RDM_TS_t; 
    end

% if size(RDM,1)==size(RDM,2)
%     RDM(logical(eye(size(RDM))))=0; % fix diagonal: zero by definition
%     RDM=squareform(RDM);
% end