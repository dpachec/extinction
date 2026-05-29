%%
%%
clear, clc
cd F:\Daniel\paris_data_extinction\paris_data_micro\included_in_2025

myFolders = dir(); myFolders=myFolders(3:end); 
myFolders= sort({myFolders.name})';

for foldi = 1:length(myFolders)

    cd([myFolders{foldi} '/micro'])


    myFiles = dir('*ncs'); myFiles = myFiles(3:end); 
    myFiles = sort({myFiles.name})';
    tmpName = cellfun(@(x) strsplit(x, '_'), myFiles, 'un', 0); 
    tmpName = unique(cellfun(@(x) x{4}, tmpName, 'un', 0)); 
    microChanNames{foldi,: } = tmpName; 

    cd ..
    cd ..


end




















%%


















%%











