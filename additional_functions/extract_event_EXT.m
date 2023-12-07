function [Ev3] = extract_event_EXT(EEG)

Ev = [{EEG.event.type}]'; 
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
clen = cellfun(@length, Ev1); 
EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
Ev2 = cat(1, Ev1{:});
Ev2(:, 10) = erase(Ev2(:, 10), ' ');
Ev3 = double(string(Ev2));


end