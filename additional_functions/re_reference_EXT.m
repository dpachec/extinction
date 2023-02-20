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
        
 
        count = 1;
        for chani = 1:length(EEG.chanlocs)
            if ~isempty(EEG.chanlocs(chani).fsLabelsR) 
                labelR = EEG.chanlocs(chani).fsLabelsR; 
                labelR = strsplit(labelR, ';');
                if sum(contains(labelR, 'White-Matter')) == 2 ...
                    | sum(contains(labelR, 'Unknown')) == 2 ...
                    | sum(contains(labelR, 'Unknown')) + sum(contains(labelR, 'White-Matter')) == 2 
                    chans2remove(count,:) = chani; 
                    count = count+ 1; 
                end
            else
                chans2remove(count,:) = chani; 
                count = count+ 1; 
            end
        end

        EEG.chanlocs(chans2remove) = []; 
        EEG.data(chans2remove, :) =  []; 

        EEG = rem_EEGLAB_fields(EEG);


       disp (' >>>> bipolar reference all electrodes');







        
    end
    
% average reference
if strcmp (mode, 'aver');
    disp ('>> AVERAGE REFERENCE...');
    dataRef = EEG.data; %dataRef(chans2exc1, :) = []; EEG.chanlocs(chans2exc1, :) = [];
    EEG_average = mean(dataRef, 1);
    EEG.data = dataRef - EEG_average;
end


if strcmp (mode, 'white_matter') %  % % CONTINUE HERE

    % % % % this loop below detects the closest white-matter channel for each channel (check EEG.chanlocs dist2WM field, to check everything is correct)
    for chani = 1:size(EEG.chanlocs,2)
        minDist = 100; 
        currCoord = EEG.chanlocs(chani).mniCoord;
        for chani2 = 1:size(EEG.chanlocs,2)
            coord2ev = EEG.chanlocs(chani2).mniCoord;
            if ~isempty(currCoord) & ~isempty(coord2ev) 
                if ~contains(EEG.chanlocs(chani).fsLabel, 'White') & contains(EEG.chanlocs(chani2).fsLabel, 'White')
                    dist2ev = norm(currCoord - coord2ev);
                    if dist2ev < minDist 
                        minDist = dist2ev;
                        id = chani2;
                        EEG.chanlocs(chani).dist2WM = minDist;
                        EEG.chanlocs(chani).idWM = id;
                        EEG.data(chani, :) = EEG.data(chani, :) - EEG.data(id, :);
                    end
                end
            end
        end
    end

    count = 1;
    for chani = 1:length(EEG.chanlocs)
        if ~isempty(EEG.chanlocs(chani).fsLabel) 
            labelR = EEG.chanlocs(chani).fsLabel; 
            if contains(labelR, 'White-Matter') | contains(labelR, 'Unknown') 
                chans2remove(count,:) = chani; 
                count = count+ 1; 
            end
        else
            chans2remove(count,:) = chani; 
            count = count+ 1; 
        end
    end

    EEG.chanlocs(chans2remove) = []; 
    EEG.data(chans2remove, :) =  []; 

    EEG = rem_EEGLAB_fields(EEG);

    % % % % only transpose for the paris data
    if sub(1) == 'p'
        EEG.chanlocs = EEG.chanlocs'; 
    end

   disp (' >>>> white-matter reference all electrodes');



end

%%end function






















