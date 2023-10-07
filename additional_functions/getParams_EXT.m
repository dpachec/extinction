
function [cfg] = getParams_EXT(f2sav)
    f2t = strsplit(f2sav, '_'); 
    cfg_prev.tyRSA = f2t{1};
    cfg_prev.roi = f2t{2};
    freqs = strsplit(f2t{4}, '-');
    cfg_prev.LT         = f2t{3}; 
    if ~( strcmp(f2t{1}, 'TR') | strcmp(f2t{1}, 'POW') )
        cfg_prev.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    else
        cfg_prev.freqs       = [];
    end
    cfg_prev.avTFV       = double(string((f2t{5}))); 
    cfg_prev.fR          = double(string((f2t{6}))); 
    tParams         = strsplit(f2t{7}, '-');
    cfg_prev.win_width   = double(string((tParams{1})));
    cfg_prev.mf          = double(string((tParams{2})));
    cfg_prev.TG          = double(string((f2t{8}))); 
    cfg_prev.contr2sav   = strsplit(f2t{9}, '-');

    if ~(cfg_prev.fR == 1 & cfg_prev.TG == 1)
        cfg = cfg_prev; 
    else
        error('FR and TG not possible together')
        return
    end
    

end