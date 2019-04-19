function analyzeTimes(experiment)


sizes = experiment.sizes;
assert(isequal(experiment.ranSizes,sizes'),"These need to be equal!")

times = experiment.times;
cpuTimes = experiment.experimentTimes;
wallTimes = experiment.wallTime;
iterations = experiment.iterations;

numSizes = length(sizes);
maxSize = max(sizes);
maxLength = size(times{1},1);
data = nan(maxLength,numSizes);
numExperiments = experiment.numExperiment;

cubic = @(x) x.^3;

figure, hold on
for s = 1:numSizes

    sAtStep = sizes(s);
    cpuTime = cpuTimes{s};
    wallTime = wallTimes{s};
    
    %cpuTies
    for e = 1:numExperiments
        xs = repmat(sAtStep,1,size(cpuTime,1));
        plot(xs,cpuTime(:,e),'ro')
        plot(xs,wallTime(:,e),'bo')
    end
    
    %Times
%     padSize = maxLength - length(timeAtSize);
%     data(:,s) = [timeAtSize;zeros(padSize,1)];
%     
% %     normedTime = timeAtSize./maxSize;
% %     meanVal = mean(normedTime);
% %     xs = repmat(log(sAtStep),1,length(timeAtSize));
% %     ys = log(timeAtSize);
%     xs = repmat(sAtStep,1,length(timeAtSize));
%     ys = timeAtSize;
% %     ys(ys == -Inf) = 0;
% % %     cubed = log(cubic(linspace(1:10000)));
%     loglog(xs,ys,'bo');
% %     plot(logspace(1:10000),cubed,'go');
end
factor = 2e-10;
plot(1:10000,cubic(1:10000).*factor,'g-');
xlabel('Size'),ylabel('Time (Seconds)');
title('Experiment Times')
legend('Cpu Time','Wall Time','x^3')
hold off

% bar graph matlab algorithm
vectors = cellfun(@transpose,cellfun(@sum,cpuTimes,'un',0),'un',0);
cpu_means = cellfun(@mean,vectors)./8;

vectors = cellfun(@transpose,cellfun(@sum,wallTimes,'un',0),'un',0);
wall_means = cellfun(@mean,vectors);

data = [cpu_means;wall_means];
figure,bar(data)

end

