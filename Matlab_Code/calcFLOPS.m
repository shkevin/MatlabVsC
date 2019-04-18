function experiment = calcFLOPS(toPlot)

if ~exist('toPlot','var') || isempty(toPlot)
    toPlot = false;
end

numExperiment = 20;
sizes = [31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,...
319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769:100:10000,10000];
iterations = flip(floor(sizes/20));
seconds = 0;

a = -1;
b = 1;

% Randomly shuffle sizes
shuffle = @(v)v(randperm(numel(v)));
randSizes = shuffle(sizes);

experiment = struct();
experiment.times = cell(length(sizes),1);
experiment.numExperiment = numExperiment;
experiment.experimentTimes = cell(length(sizes),1);
experiment.wallTime = cell(length(sizes),1);
experiment.sizes = sizes;
experiment.iterations = iterations;
experiment.runTime = seconds;
experiment.ranSizes = nan(length(sizes),1);

if toPlot
    figure, hold on
end

total = tic;
for s = randSizes
    seconds = 0;
    sizeIndex = find(sizes == s);
    
    %Generate the matrix A and matrix B outside cputime.
    A = (b-a).*rand(s) + a;
    B = (b-a).*rand(s) + a;

    % Warm up the cache first
    A*B;
    
    iterationSize = iterations(sizeIndex);
    times = zeros(iterationSize,1);
    wallTimes = zeros(iterationSize,numExperiment);
    experimentTimes = zeros(iterationSize,numExperiment);
    
    for i = 1:iterationSize
  
        % Compute the solution for A*B numExperiment times.
        for j = 1:numExperiment
            tic;
            expTime = cputime;
            A*B;
            elapsed = cputime-expTime;
            walltime = toc;
            wallTimes(i,j) = walltime;
            experimentTimes(i,j) = elapsed;
        end
%         iterationTime = (cputime-iterTime);
        iterationTime = sum(experimentTimes(i,:));
        seconds = seconds + iterationTime;
        times(i) = iterationTime;
        
        if toPlot
            plot(s,iterationTime,'k-*','LineWidth',1);
        end
    end
    
    experiment.wallTime{sizeIndex} = wallTimes;
    experiment.experimentTimes{sizeIndex} = experimentTimes;
    experiment.times{sizeIndex} = times;
    experiment.ranSizes(sizeIndex) = s;
end

experiment.runTime = toc(total);
% x = [1 iterations];
% p = polyfit(x,times',4);
% f1 = polyval(p,x);
% plot(x,f1,'r--','linewidth',1.5)
if toPlot
    xlabel('Input Size');ylabel('Time (Seconds)');
    title('Asymptotic Growth')
    hold off
end

%Save the simulation based off aquisition date.
defaultName = datetime;
defaultName = strtrim(regexprep(string(defaultName),{'-',' ',':'},{'','_',''}));

defaultFolder = '/users/shkevin/assignment3/Matlab_Code';
fullfileName = strcat(defaultFolder,defaultName);
save(strcat(fullfileName,'.mat'), 'experiment');

end

