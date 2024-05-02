
function [cfg] = getParams_EXT(f2sav)
    f2t = strsplit(f2sav, '_'); 
    

    switch f2t{1}
        case 'nRSM'
            cfg.roi = f2t{2};
            cfg.LT         = f2t{3}(1); 
            cfg.period     = f2t{3}(2);
            freqs = strsplit(f2t{4}, '-');
            cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
            
            cfg.avTFV       = double(string((f2t{5}))); 
            cfg.fR          = double(string((f2t{6}))); 
            tParams         = strsplit(f2t{7}, '-');
            cfg.win_width   = double(string((tParams{1})));
            cfg.mf          = double(string((tParams{2})));
            cfg.powF2load   = ['POW_' cfg.roi '_' cfg.LT '_' f2t{4} '_' num2str(cfg.mf)]; 

        case 'POW'
            cfg.roi         = f2t{2};
            cfg.LT          = f2t{3}; 
            freqs = strsplit(f2t{4}, '-');
            cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
            cfg.tR          = double(string((f2t{5})));
    
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