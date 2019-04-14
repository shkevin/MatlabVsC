function experiment = calcFLOPS(toPlot)

if ~exist('toPlot','var') || isempty(toPlot)
    toPlot = false;
end

iterations = 1:10;
numExperiment = 100;
sizes = [31 32 96 97 127 128 129 191 192 229 255 256 257 319 320 321 417 ...
    479 480 511 512 639 640 767 768 769];
seconds = 0;

a = -1;
b = 1;

% Randomly shuffle sizes
shuffle = @(v)v(randperm(numel(v)));
sizes = shuffle(sizes);

experiment = struct();
experiment.times = nan(length(iterations),length(sizes));
experiment.sizes = sizes;
experiment.runTime = seconds;

if toPlot
    figure, hold on
end

sizeIndex = 1;
tic;
for s = sizes
    times = zeros(length(iterations),1);
    seconds = 0;
    
    %Generate the matrix A and matrix B outside cputime.
    A = (b-a).*rand(s) + a;
    B = (b-a).*rand(s) + a;

    % Warm up the cache first
    A*B;
    
    for i = iterations
        
        t = cputime;
        
        % Compute the solution for A,b 10 times.
        for j = 1:numExperiment
            A*B;
        end
        iterationTime = (cputime-t);
        seconds = seconds + iterationTime;
        averageTime = iterationTime/numExperiment;
        times(i) = averageTime;
        
        if toPlot
            plot(s,averageTime,'k-*','LineWidth',1);
        end
    end
    experiment.times(:,sizeIndex) = times;
    sizeIndex = sizeIndex + 1;
end

experiment.runTime = toc;
% x = [1 iterations];
% p = polyfit(x,times',4);
% f1 = polyval(p,x);
% plot(x,f1,'r--','linewidth',1.5)
if toPlot
    xlabel('Input Size');ylabel('Time (Seconds)');
    title('Asymptotic Growth')
    hold off
end

end

