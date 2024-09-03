
function [cfg] = getParams_EXT(f2sav)
    f2t = strsplit(f2sav, '_'); 
    

    switch f2t{1}
        case 'nRSM'
            cfg.roi = f2t{2};
            cfg.LT         = f2t{3}(1); 
            cfg.period     = f2t{3}(2:end);
            freqs = strsplit(f2t{4}, '-');
            cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
            
            cfg.avTFV       = double(string((f2t{5}))); 
            cfg.fR          = double(string((f2t{6}))); 
            tParams         = strsplit(f2t{7}, '-');
            cfg.win_width   = double(string((tParams{1})));
            cfg.mf          = double(string((tParams{2})));
            %cfg.powF2load   = ['POW_' cfg.roi '_' cfg.LT '_' num2str(cfg.mf/10)]; 
            cfg.powF2load   = ['POW_' cfg.roi '_' cfg.LT '_10']; 

        case 'POW'
            cfg.roi         = f2t{2};
            cfg.LT          = f2t{3}; 
            cfg.tR          = double(string((f2t{4})));
        case 'RSA'
            cfg.roi = f2t{2};
            cfg.LT         = f2t{3}; 
            freqs = strsplit(f2t{4}, '-');
            cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
            
            cfg.avTFV       = double(string((f2t{5}))); 
            cfg.fR          = double(string((f2t{6}))); 
            tParams         = strsplit(f2t{7}, '-');
            cfg.win_width   = double(string((tParams{1})));
            cfg.mf          = double(string((tParams{2})));
            %cfg.powF2load   = ['POW_' cfg.roi '_' cfg.LT '_' num2str(cfg.mf/10)]; 
            cfg.powF2load   = ['POW_' cfg.roi '_' cfg.LT '_10']; 
            cfg.TG = double(string(f2t{8})); 
            cfg.contr2sav   = strsplit(f2t{9}, '-');


    
end


% 
% function [cfg] = getParams_EXT(f2sav)
%     f2t = strsplit(f2sav, '_'); 
%     cfg.tyRSA = f2t{1};
%     cfg.roi = f2t{2};
%     freqs = strsplit(f2t{4}, '-');
%     cfg.LT         = f2t{3}; 
%     if ~( strcmp(f2t{1}, 'TR') )
%         cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
%     else
%         cfg.freqs       = [];
%     end
%     cfg.avTFV       = double(string((f2t{5}))); 
%     cfg.fR          = double(string((f2t{6}))); 
%     tParams         = strsplit(f2t{7}, '-');
%     cfg.win_width   = double(string((tParams{1})));
%     cfg.mf          = double(string((tParams{2})));
%     cfg.TG          = double(string((f2t{8}))); 
%     cfg.contr2sav   = strsplit(f2t{9}, '-');
% 
% 
%     if strcmp(f2t{1}, 'TR') & cfg.fR ==1
%         error('traces contrast cannot be Frequency resolved')
%     end
% 
%     if cfg.fR == 1 & cfg.TG == 1
%         error('FR and TG not possible together')
%     end
% 
% 
% 
% end