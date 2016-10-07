%Scanning Script
clear
close all
clc

newClouds = ScanAndRegisterFromFiles('kuukpak scan models.txt','kuukpak planar model.txt','kuukpak sphere models.txt',1);

%% Plot Clouds in common coodinate system
if 0
    figure, hold on
    cols = rand(length(newClouds),3);
    for i = 1:length(newClouds)
        newClouds(i).plotCloud3D('reg','old',cols(i,:));
    end
    axis equal
    title('Registered Coordinates')
end
%%
for i = 1:length(newClouds)
    disp(['Cloud number: ', num2str(i)])
    disp('Registration parameters in meters and degrees: (x, y, z, omega, phi, kappa)');
    disp(newClouds(i).scan.regParams(1:3)');
    disp(newClouds(i).scan.regParams(4:6)'*180/pi);
    disp('Registration parameter standard deviations: (mm, arc-seconds)')
    disp(sqrt(diag(newClouds(i).scan.regParamVarCovar(1:3,1:3)))'*1000);
    disp(sqrt(diag(newClouds(i).scan.regParamVarCovar(4:6,4:6)))'*3600*180/pi);
    disp(' ')
end

%%
allVariances = [];
RegErrs = [];

for i = 1:length(newClouds)
    allVariances = [allVariances; newClouds(i).varRegXYZ];
    disp(['Mean Standard Deviation of X, Y, Z post registration (mm) for cloud number ', num2str(i), ':'])
    disp(mean(sqrt(newClouds(i).varRegXYZ))*1000)
    RegErrs = [RegErrs;  newClouds(i).maxRegSTD];
end

disp('Mean Standard Deviation of X, Y, Z post registration for all clouds(mm):')
disp(mean(sqrt(allVariances))*1000)


%%
figure, subplot(1,2,1), hist(RegErrs,100), title('Max Errors per point'), xlabel('Meters')
subplot(1,2,2), hist(log10(RegErrs),100), title('Log of Max Errors per point'), xlabel('Log of Meters')

%%

% d = displayCloudsFromArray(newClouds,10);

%% output pts

if 1
    fid = fopen('TLSPointCloud.pts','w');
    
    for i =  1:length(newClouds)
        for j = 1:length(newClouds(i).GLOBALXYZ)
            outputString = [num2str(newClouds(i).GLOBALXYZ(j,1)),', ',num2str(newClouds(i).GLOBALXYZ(j,2)),...
                ', ',num2str(newClouds(i).GLOBALXYZ(j,3)),', ',num2str(newClouds(i).maxRegSTD(j)),'\n'];
            if newClouds(i).maxRegSTD(j) ~= 0;
                fprintf(fid,outputString);
            end
        end
    end
    fclose(fid);
end