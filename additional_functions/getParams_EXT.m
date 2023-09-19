
function [cfg] = getParams_EXT(f2sav)
    f2t = strsplit(f2sav, '_'); 
    freqs = strsplit(f2t{1}, '-');
    cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    cfg.avTFV       = double(string((f2t{2}))); 
    cfg.fR          = double(string((f2t{3}))); 
    cfg.fitMode     = double(string((f2t{4}))); %1= trials; 0 = no trials; 
    tParams         = strsplit(f2t{5}, '-');
    cfg.win_width   = double(string((tParams{1})));
    cfg.mf          = double(string((tParams{2})));
    cfg.TG          = string((f2t{6})); 
    cfg.contr2sav   = strsplit(f2t{7}, '-');
    

end