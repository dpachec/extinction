function EEG = re_reference_WM(EEG, mode, sub)

disp ('>>>>> re-referencing...');

%first remove EEG and POL for chinese data
chs = {EEG.chanlocs(:).labels}';
chs = erase(chs, 'POL'); %deletes over the whole array
chs = erase(chs, 'EEG');
chs = erase(chs, '-Ref');
chs = erase(chs, 'G_');
chs = erase(chs, '-');
chs = erase(chs, ' ');
chs = erase(chs, '''');


if strcmp (mode, 'bipo');

        disp ('>> BIPOLAR REFERENCE...');

        if sub(1) == 'p'
            chp1 = cellfun(@(x) strsplit(x, '_'), chs, 'un', 0);
            chp1 = cellfun(@(x) x{1}, chp1, 'un', 0);
            chp1 = strrep(chp1, '1', 'ONE'); 
            chp1 = strrep(chp1, '2', 'TWO');
            chp1 = strrep(chp1, '3', 'THREE');
            chp1 = strrep(chp1, '4', 'FOUR');
            chp1 = strrep(chp1, '5', 'FIVE');
            chp1 = strrep(chp1, '6', 'SIX');
            chp1 = strrep(chp1, '7', 'SEVEN');
            chp1 = strrep(chp1, '8', 'EIGHT');
            chp1 = strrep(chp1, '9', 'NINE');
            chp1 = strrep(chp1, '0', 'ZERO');
            chp2 = cellfun(@(x) strsplit(x, '_'), chs, 'un', 0);
            clen = cellfun(@length, chp2); 
            chp2(clen==1)= strcat(chp2(clen==1), '_X');
            
            chp2 = cellfun(@(x) x{2}, chp2, 'un', 0);
            chs = strcat(chp1,chp2);
            
        end

        for i = 1:length(EEG.chanlocs)-1


            currChan = chs{i};nextChan = chs{i+1};
            index = find(isletter(currChan));currChan = currChan(index);
            index = find(isletter(nextChan));nextChan = nextChan(index);
            currChanB = chs{i};nextChanB = chs{i+1};
            indexNum = find(~isletter(currChanB)); currNum = str2num(currChanB(indexNum));
            indexNumNext = find(~isletter(nextChanB)); nextNum = str2num(nextChanB(indexNumNext));
            
            disp([currChan ' ' nextChan]);
            disp(['number >> ' num2str(nextNum) ' ' num2str(currNum)]);

            if  strcmp(currChan, nextChan) & (nextNum - currNum == 1)
                disp(['Equal: ' EEG.chanlocs(i).labels ' = ' EEG.chanlocs(i+1).labels ' - ' EEG.chanlocs(i).labels]);
                EEG.data(i,:) = EEG.data(i+1,:,:)  - EEG.data(i,:,:);
                if ~isempty(EEG.chanlocs(i).mniCoord) & ~isempty(EEG.chanlocs(i+1).mniCoord) 
                    mniCoordR = (EEG.chanlocs(i).mniCoord + EEG.chanlocs(i+1).mniCoord) / 2; 
                    EEG.chanlocs(i).mniCoordR = mniCoordR; 
                    EEG.chanlocs(i).labelsR = [EEG.chanlocs(i).labels ';' EEG.chanlocs(i+1).labels ];
                    EEG.chanlocs(i).fsLabelsR = [EEG.chanlocs(i).fsLabel ';' EEG.chanlocs(i+1).fsLabel ];
                end
            else
                disp('Different');
                ids2rem(i) = 1;
            end 
        end
        EEG.data(end,:) = []; %remove last channel (not bipolarized)
        EEG.chanlocs(end) = [];
        ids2rem = logical (ids2rem);
        EEG.chanlocs(ids2rem) = []; 
        EEG.data(ids2rem, :) = [];
        
        disp (' >>>> bipolar reference all electrodes');














        
    end
    
% average reference
if strcmp (mode, 'aver');
    disp ('>> AVERAGE REFERENCE...');
    dataRef = EEG.data; %dataRef(chans2exc1, :) = []; EEG.chanlocs(chans2exc1, :) = [];
    EEG_average = mean(dataRef, 1);
    EEG.data = dataRef - EEG_average;
end

%%end function






















