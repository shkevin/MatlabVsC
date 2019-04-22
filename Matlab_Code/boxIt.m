function boxIt(data)
% Function:
%   boxIt(data):
%       Produces a box plot of the given data set.
%
% References:
%	https://www.sciencedirect.com/topics/medicine-and-dentistry/rank-sum-test
%   https://www.mathworks.com/matlabcentral/answers/426528-indicating-statistical-significance-on-boxplot-in-matlab       
%
% Input:
%   data: Matrix of row vectors of data of interest.
%
% Output:
%   -None-


figure,
boxplot(data,'Labels',{'Matlab','genvect'},'notch','on');
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
axis([xlim,4.5,ceil(max(yt)*1.2)])

hold on
% Plot the line
plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k')

% Plot the stars
center = mean(xt([1 2]));
markers = [mean(xt([1 2]))-0.05,center,mean(xt([1 2]))+0.05,mean(xt([1 2]))+0.1];
plot(markers ,max(yt)*1.15, '*k')

title('Matlab-Genvect Boxplot')
xlabel('n = 3000'),ylabel('Time (seconds)');
hold off

end

