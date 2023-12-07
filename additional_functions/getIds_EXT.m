function [ids] = getIds_EXT(Ev2, cfg)

    Ev3 = double(string(Ev2)); 
    if strcmp(cfg.contr2sav{1}, 'ALLA')
        ids = Ev3(:, 2) == 1;
    end
    if strcmp(cfg.contr2sav{1}, 'ALLE')
        ids = Ev3(:, 2) == 2;
    end

    if strcmp(cfg.contr2sav{1}, 'CSPA')
        ids = Ev3(:, 2) == 1 & Ev3(:, 8) == 1 ;
    end
    if strcmp(cfg.contr2sav{1}, 'CSMA')
        ids = Ev3(:, 2) == 1 & Ev3(:, 8) == 0;
    end
    if strcmp(cfg.contr2sav{1}, 'CSPE')
        ids = Ev3(:, 2) == 2 & Ev3(:, 8) == 1 ;
    end
    if strcmp(cfg.contr2sav{1}, 'CSME')
        ids = Ev3(:, 2) == 2 & Ev3(:, 8) == 0;
    end





end