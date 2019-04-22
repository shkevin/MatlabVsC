function wilcoxonRankSum(data)
% Function:
%   wilcoxonRankSum(data):
%       Performs a Wilcoxon Rank-Sum Test on the data.
%
% References:
%       https://www.sciencedirect.com/topics/medicine-and-dentistry/rank-sum-test
%
% Input:
%   data: Matrix of row vectors of data of interest.
%
% Output:
%   -None-

%define all data in table
assert(size(data,2) > 1,"Can't compare with only one data set");
data = data;

% Assumed square data set since table;
n = size(data,1);

% Calculate ranks
calcMean = @(n1,n2) n1*(n1+n2+1)/2;
calcStd = @(n1,n2) sqrt(n1*n2*(n1+n2+1)/n);
ranks = rankIt(data);

totalSum = sum(ranks(:,1));
meanVal = calcMean(n,n);
stdVal = calcStd(n,n);

[p,h,stat] = ranksum(data(:,1),data(:,2))
% testStat = (totalSum - meanVal)/stdVal;

% assert(stat.ranksum == totalSum,"Sums need to be equal");
testStat = stat.zval;

n = size(data(:,1),1);
lower = round((n/2) - 1.96*(sqrt(n)/2));
upper = round(1+(n/2)+1.96*(sqrt(n)/2));
testData = data(:,2);
% Always assumes column 1 is the data to compare against
confidenceInterval = [testData(lower),testData(upper)]

if h == 0
    disp('Fail to reject the null hypothesis - The medians are statistically equal!')
else
    disp('Reject the null hypothesis - The median values are not equal!')
end

end

