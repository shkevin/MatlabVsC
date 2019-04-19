function [cputimes,walltimes] = readCData(algorithm)

base = @(times)strcat('./CS533_Assignment3_Matmult/',algorithm,'/experimentInfo/',times);
cpuPath = base('cputime');
wallPath = base('walltime');

% Assumes all folders are the same.
cpu_files = dir(fullfile(cpuPath,'*.txt'));
wall_files = dir(fullfile(wallPath,'*.txt'));

assert(isequal(length(cpu_files),length(wall_files)));

n = length(cpu_files);
data = cell(n,4);
cpu_sizes = nan(n,1);
wall_sizes = nan(n,1);
for i = 1:n
    cpu_file = cpu_files(i).name;
    wall_file = wall_files(i).name;
    
    
    cputimeAtI = load(fullfile(cpuPath,cpu_file));
    walltimeAtI = load(fullfile(wallPath,wall_file));
    
    data{i,1} = cputimeAtI;
    data{i,3} = walltimeAtI;
    data{i,2} = str2double(regexprep(cpu_files(i).name,'.txt',''));
    data{i,4} = str2double(regexprep(wall_files(i).name,'.txt',''));
end

cputimes = sort(data{:,1:2},2);
walltimes = sort(data{:,3:4},2);

end

