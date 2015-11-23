function [h1,h2] = displayCloudsFromArray(cloudArray)

regXYZ = [];
varRegXYZ = [];
for i = 1:length(cloudArray)
    regXYZ = [regXYZ; cloudArray(i).regXYZ];
    varRegXYZ = [varRegXYZ; cloudArray(i).varRegXYZ];    
end

totalDev = sqrt(sum(varRegXYZ,2));

data = sortrows([totalDev,regXYZ]);

numOfCategories = 255;

indexes = round( 1:length(data)/(numOfCategories):length(data))';

%% Histogram Equalization
EqualizedColours = (indexes-1)/max(indexes-1);
equalizedCols = sqrt([EqualizedColours ,zeros(size(EqualizedColours)),1-EqualizedColours]);

h1 = figure;
hold on

for i = 2:numOfCategories
    plot3(data(indexes(i-1):indexes(i),2),data(indexes(i-1):indexes(i),3),...
        data(indexes(i-1):indexes(i),4),'.','Color',equalizedCols(i,:));
end
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Equalized Histogram')

%% Normalized colours
normColours = data(:,1)-data(1,1);
normColours = data(:,1)./max(normColours);
normCol = sqrt([normColours ,zeros(size(normColours)),1-normColours]);
h2 = figure;
hold on

for i = 2:numOfCategories
    plot3(data(indexes(i-1):indexes(i),2),data(indexes(i-1):indexes(i),3),...
        data(indexes(i-1):indexes(i),4),'.','Color',normCol(indexes(i-1),:))
end
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Normal Histogram')
