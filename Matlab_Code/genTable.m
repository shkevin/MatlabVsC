function [T,best,speedups] = genTable()

format short
cols = {'Size','Matlab','Row','Blas','Col','Genvect','Autovect','Naive','Copy','Blocked'};

sizes = [31,32,96,97,127,128,129,191,192,229,255,256,257,319,320,321,417,...
    479,480,511,512,639,640,767,768,769,1000,3000,4000,10000];
tableSize = [length(sizes) length(cols)];

types = {'double','double','double','double','double','double','double','double','double','double'};
T = table('Size',tableSize,'VariableTypes',types,'variableNames',cols);

T{:,1} = sizes';

base = './Matlab_Code/AllData/';

matlab = load(fullfile(base,'matlab'));
row = load(fullfile(base,'row'));
blas = load(fullfile(base,'blas'));
col = load(fullfile(base,'col'));
genvect = load(fullfile(base,'genvect'));
autovect = load(fullfile(base,'autovect'));
naive = load(fullfile(base,'naive'));
copy = load(fullfile(base,'copy'));
blocked = load(fullfile(base,'blocked'));

meanIt = @(data) cellfun(@mean,cellfun(@mean, data,'un',0));

T.Matlab = meanIt(matlab.matlab.experimentTimes);
T.Row = meanIt(row.row.experimentTimes);
T.Blas = meanIt(blas.blas.experimentTimes);
T.Col = meanIt(col.col.experimentTimes);
T.Genvect = meanIt(genvect.genvect.experimentTimes);
T.Autovect = meanIt(autovect.autovect.experimentTimes);
T.Naive = meanIt(naive.naive.experimentTimes);
T.Copy = meanIt(copy.copy.experimentTimes);
T.Blocked = meanIt(blocked.blocked.experimentTimes);

m = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
% colors = {[0, 0, 1],[1, 0, 0],[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0, 0.5, 0],[0, 0.4470, 0.7410],[0.75, 0, 0.75],[0.4940, 0.1840, 0.5560],[0.3010, 0.7450, 0.9330]};
xs = sizes;
cubic = @(x) x.^3;
figure,hold on
for alg = 2:length(cols)
    ys = T{:,alg};
    h = plot(xs,ys,...
        m{alg-1},...
        'MarkerSize',10,...
        'LineWidth',1);
end
factor = @(x) mean(x{end})/(max(sizes)^3);

colors = get(gca,'colororder');
xs = 1:10000;
ys = cubic(1:10000);
factors = cell(9,1);
factors{1} = factor(matlab.matlab.experimentTimes);
factors{2} = factor(row.row.experimentTimes);
factors{3} = factor(blas.blas.experimentTimes);
factors{4} = factor(col.col.experimentTimes);
factors{5} = factor(genvect.genvect.experimentTimes);
factors{6} = factor(autovect.autovect.experimentTimes);
factors{7} = factor(naive.naive.experimentTimes);
factors{8} = factor(copy.copy.experimentTimes);
factors{9} = factor(blocked.blocked.experimentTimes);

plot(xs,ys.*factors{1},'Color',colors(2,:))
plot(xs,ys.*factors{2},'Color',colors(2,:))
plot(xs,ys.*factors{3},'Color',colors(3,:))
plot(xs,ys.*factors{4},'Color',colors(4,:))
plot(xs,ys.*factors{5},'Color',colors(5,:))
plot(xs,ys.*factors{6},'Color',colors(6,:))
plot(xs,ys.*factors{7},'Color',colors(7,:))
plot(xs,ys.*factors{8},'Color',colors(1,:))
plot(xs,ys.*factors{9},'Color',colors(2,:))

title('Asymptotic Growth of all Algorithms')
xlabel('Size'),ylabel('Time (seconds)');
legend(cols(2:end),'Location','west');
hold off

[~,idx] = min(T{:,:},[],2);
bestAlgs = unique(idx);
sampleIdx = [1:5:30, 27:30];

best = table();
best.Size = T.Size(sampleIdx);
best.Matlab = T.Matlab(sampleIdx,:);
best.Copy = T.Copy(sampleIdx,:);
best.Genvect = T.Genvect(sampleIdx,:);

speedups = table();
speedups.Algorithm = cols(2:end)';
speedups.Growth_Factor = factors;
speedups.Speed_Up = cell2mat(factors(:))./factors{1};

% vectors = @(data) cellfun(@mean, cellfun(@transpose,cellfun(@sum,data,'un',0),'un',0));
meanIt = @(x) cellfun(@mean,cellfun(@mean, x,'un',0));
padthedata = @(x) [x;nan(2500-size(x,1),1)];

index = 1;
for i = sampleIdx
    vals = matlab.matlab.experimentTimes{i};
    data(:,index) = padthedata(vals(:));
    index = index + 1;
end

boxplot(data,'Labels',sizes(sampleIdx))

end

