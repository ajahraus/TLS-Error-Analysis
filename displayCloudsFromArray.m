function indexes = displayCloudsFromArray(cloudArray,numOfCategories)

regXYZ = [];
maxRegSTD = [];
for i = 1:length(cloudArray)
    regXYZ = [regXYZ; cloudArray(i).regXYZ];
    maxRegSTD = [maxRegSTD; cloudArray(i).maxRegSTD];    
end

data1 = sortrows([maxRegSTD,regXYZ]);

normalizedDistance = abs((data1(:,1) - mean(data1(:,1)))./std(data1(:,1)));
key = normalizedDistance<3;
data = deleteRowKey(data1,key);
outliers = deleteRowKey(data1,~key);


% disp(data(1,1))
% disp(data(end,1))

indexes = 1:length(data)/(numOfCategories-1):length(data);
indexes = round([indexes, length(data)])';

% Histogram Equalization
EqualizedColours = 0.67 -  (indexes-1)/(max(indexes-1))*0.67;
temp = ones(length(EqualizedColours),1,3);
temp(:,1,1) = EqualizedColours;
rgb = hsv2rgb(temp);
equalizedCols(:,3) = rgb(:,1,3);
equalizedCols(:,2) = rgb(:,1,2);
equalizedCols(:,1) = rgb(:,1,1);
% equalizedCols = sqrt([EqualizedColours ,zeros(size(EqualizedColours)),1-EqualizedColours]);

h1 = figure;
hold on


legendStrings = cell(numOfCategories,1);
numOdigits = 3;
for i = 2:numOfCategories-1
    figure(h1);
    plot3(data(indexes(i-1):indexes(i),2),data(indexes(i-1):indexes(i),3),...
        data(indexes(i-1):indexes(i),4),'.','Color',equalizedCols(i-1,:));
%     plot(data(indexes(i-1):indexes(i),2),data(indexes(i-1):indexes(i),3),...
%         '.','Color',equalizedCols(i,:));

    legendStrings{i-1} = [num2str(data(indexes(i-1),1)*1000,numOdigits),' : ',num2str(1000*data(indexes(i),1),numOdigits)];
    
%     figure(h2)
%     histogram(data(indexes(i-1):indexes(i),1),50,'Color',equalizedCols(i-1,:));
    
    if i == numOfCategories-1
        figure(h1)
        plot3(data(indexes(i):indexes(end),2),data(indexes(i):indexes(end),3),...
        data(indexes(i):indexes(end),4),'.','Color',equalizedCols(i,:))
        legendStrings{i} = [num2str(1000*data(indexes(i),1),numOdigits), ' : ',num2str(1000*data(indexes(end),1),numOdigits)];
        
%         figure(h2)
%         histogram(data(indexes(i):indexes(end),1),50,'Color',equalizedCols(i,:));
    end
    
end

plot3(outliers(:,2),outliers(:,3),outliers(:,4),'.','Color',[1,0,0])
legendStrings{end} = [num2str(1000*min(outliers(:,1)),numOdigits),' : ',num2str(1000*max(outliers(:,1)),numOdigits)];

figure(h1)
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Equalized Histogram')
legend(legendStrings,'Location','EastOutside')


% h2 = figure;
% hold on
% figure(h2)
% plot3(outliers(:,2),outliers(:,3),outliers(:,4),'r.')
% axis equal
% title('Outliers')
