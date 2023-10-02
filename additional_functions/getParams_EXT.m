
function [cfg] = getParams_EXT(f2sav)
    f2t = strsplit(f2sav, '_'); 
    freqs = strsplit(f2t{2}, '-');
    if ~strcmp(f2t{1}, 'T') 
        cfg.freqs       = [double(string((freqs{1}))) : double(string((freqs{2})))];
    else
        cfg.freqs       = [];
    end
    cfg.avTFV       = double(string((f2t{3}))); 
    cfg.fR          = double(string((f2t{4}))); 
    cfg.fitMode     = double(string((f2t{5}))); %1= trials; 0 = no trials; 
    tParams         = strsplit(f2t{6}, '-');
    cfg.win_width   = double(string((tParams{1})));
    cfg.mf          = double(string((tParams{2})));
    cfg.TG          = double(string((f2t{7}))); 
    cfg.contr2sav   = strsplit(f2t{8}, '-');
    

end