function [factor,difference] = testIt(experiment)
% Function:
%   [factor,difference] = testIt(experiment):
%       More of a throw-away function that plots the means of the given
%       experiment. This was just used to help understand underlying data.   
%
% Input:
%   experiment: Struct containing data and details from experiment.
%
% Output:
%   factor: Growth factor of the given experiment.
%   difference: Differences between wall time and CPU time.

cputimes = experiment.experimentTimes;
walltimes = experiment.wallTimes;
sizes = experiment.sizes;
numSizes = length(sizes);
cubic = @(x) x.^3;

% for s = 1:numSizes
%     cputime = cputimes{s};
%     walltime = walltimes{s};
%     
%     xs = repmat(sizes(s),1,size(cputime,1));
%     plot(xs,cputime,'rs')
%     plot(xs,walltime,'bo')
% end
% legend('Cpu Time','Wall Time','Location','west')
factor = mean(cputimes{end})/(max(sizes)^3);
% plot(1:10000,cubic(1:10000).*factor,'g-','DisplayName',sprintf('cubic reference - factor %d',factor));

vars = cellfun(@var,cellfun(@var, cputimes,'un',0),'un',0);
totalVariances = sum(cell2mat(vars))

figure;hold on

meanWall = cellfun(@mean,cellfun(@mean, walltimes,'un',0));
meanCpu = cellfun(@mean,cellfun(@mean,cputimes,'un',0));
plot(sizes,meanCpu,'rs')
plot(sizes,meanWall,'bo');
legend('CPU Times','Wall Times')
plot(1:10000,cubic(1:10000).*factor,'--r','DisplayName','Cubic reference')
plot(1:10000,cubic(1:10000).*factor./8,'-b','DisplayName','Cubic reference - divided 8')
difference = meanWall-meanCpu;
legend('Location','west')

title('Asymptotic Growth')
xlabel('Size'),ylabel('Time (seconds)')
hold off

end

