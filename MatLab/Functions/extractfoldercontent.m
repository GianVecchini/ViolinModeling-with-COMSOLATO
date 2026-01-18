function [dirContent, n] = extractfoldercontent(path)

dirListing = dir(path);         %list directory content
dirContent = {dirListing.name}; % getting data content names

if length(dirContent) <= 2
    error('Data directory is empty or does not exist');
end

dirContent = dirContent(3:end); %discard '.' and '..' directory references 
n = length(dirContent);

end