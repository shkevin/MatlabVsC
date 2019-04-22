function ranks = rankIt(data)
% Function:
%   ranks = rankIt(data):
%       Calculate the ranks for a Wilcoxon Rank-Sum test.
%
% References:
%   https://www.mathworks.com/matlabcentral/answers/165663-how-to-give-ranking-from-highest-to-lowest [1]
%   https://www.mathworks.com/matlabcentral/answers/336500-finding-the-indices-of-duplicate-values-in-one-array [2]    
%
% Input:
%   data: data to rank.
%
% Output:
%   ranks : ranks of the given data set.

% Pulled from [1]
dims = size(data);
data = data(:);
[~,p] = sort(data);
r = 1:length(data);
r(p) = r;

% Handle repeated values
% Pulled from: [2]
[~, uniqueIdx] = unique(data);
duplicates = ismember(data,find(data(setdiff(1:numel(data), uniqueIdx))));

% Determine if there were midranks that need to be calculated.
if ~isempty(duplicates)
    r(duplicates) = min(r(duplicates));
    r = midRanks(r);
end
ranks = reshape(r,dims);

end

function r = midRanks(x)
% function r = midRanks(x)
% P.B. Stark, statistics.berkeley.edu/~stark  2/17/2003
% vector of midranks of a vector x.
% midrank is (#<=) - ((#=) - 1)/2 = ((#<=) + (#<) + 1)/2
for j = 1:length(x)
    r(j) = sum(x <= x(j)) - (sum(x == x(j)) - 1)/2;
end
end