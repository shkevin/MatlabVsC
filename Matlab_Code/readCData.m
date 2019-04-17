function data = readCData(dirPath)

% Assumes all folders are the same.
files = dir(fullfile(dirPath,'*.txt'));

n = length(files);
data = cell(n,2);
for i = 1:1
    file = files(i).name;
    matrix = load(fullfile(dirPath,file));
    
    data{i,1} = matrix;
    data{i,2} = regexprep(files(i).name,'.txt','');
end

data = sort(data,2);

end

