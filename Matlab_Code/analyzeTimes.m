function analyzeTimes(experiment)


sizes = experiment.sizes;
assert(isequal(experiment.ranSizes,sizes),"These need to be equal!")

times = experiment.times;
cpuTimes = experiment.experimentTimes;
wallTimes = experiment.wallTime;
iterations = experiment.iterations;

numSizes = length(sizes);
maxSize = max(sizes);
maxLength = size(times{1},1);
data = nan(maxLength,numSizes);

figure, hold on
for s = 1:numSizes
    timeAtSize = times{s};
    sAtStep = sizes(s);
    
    padSize = maxLength - length(timeAtSize);
    data(:,s) = [timeAtSize;zeros(padSize,1)];
    
%     normedTime = timeAtSize./maxSize;
%     meanVal = mean(normedTime);
%     xs = repmat(log(sAtStep),1,length(timeAtSize));
%     ys = log(timeAtSize);
    xs = repmat(sAtStep,1,length(timeAtSize));
    ys = timeAtSize;
%     ys(ys == -Inf) = 0;
% %     cubed = log(cubic(linspace(1:10000)));
    loglog(xs,ys,'bo');
%     plot(logspace(1:10000),cubed,'go');
end

hold off


end

