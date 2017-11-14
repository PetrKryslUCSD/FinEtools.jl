function M = readacsv(filename)
filename = [filename '.csv'];
delimiter = ',';

fileID = fopen(filename,'r');
M = dlmread(filename,delimiter) 
fclose(fileID);
end