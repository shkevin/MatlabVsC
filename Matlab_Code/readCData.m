function experiment = readCData(algorithm)
% Function:
%   experiment = readCData(algorithm)
%       Reads in the desired C data simulation.
%
% Input:
%   algorithm: String containing which algorithm folder to read in.
%
% Output:
%   experiment: Struct containing data and experiment details.


base = @(times)strcat('./CS533_Assignment3_Matmult/',algorithm,'/experimentInfo/',times);
cpuPath = base('cputime');
wallPath = base('walltime');

% Assumes all folders are the same.
cpu_files = dir(fullfile(cpuPath,'*.txt'));
wall_files = dir(fullfile(wallPath,'*.txt'));

% assert(isequal(length(cpu_files),length(wall_files)));

n = length(cpu_files);
data = cell(n,4);

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

data = sortrows(data,2);

% assert(isequal(data(:,2),data(:,4)),"sizes need to be the same")

experiment = struct();
experiment.experimentTimes = data(:,1);
experiment.wallTimes = data(:,3);
experiment.sizes = cell2mat(data(:,2));

end

