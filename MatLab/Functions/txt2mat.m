function Data = txt2mat(imported_filename)
% TXT2MAT
% This function import the frf data from .txt file into matlab.
   


addpath('./data')
new_filename = strcat(strrep(imported_filename,'.txt',''),'_clean.txt');

Data = fileread(imported_filename);

Data = strrep(Data, ',', '.');
fID = fopen(new_filename,'w');
fwrite(fID,Data,'char');
fclose(fID);

Data = dlmread(new_filename,'');
delete(new_filename);