function[files] = getFiles(path)

files = dir([path]); cbytes = [files.bytes]; 
files = files(cbytes > 100);files = sort({files.name})';

end