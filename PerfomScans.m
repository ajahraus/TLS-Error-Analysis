% Scanning Script
clear
close all
clc

newClouds = ScanAndRegisterFromFiles('kuukpak scan models.txt','kuukpak planar model.txt','kuukpak sphere models.txt');


%% Plot Clouds in common coodinate system
figure, hold on
cols = rand(length(newClouds),3);
for i = 1:length(newClouds)
    newClouds(i).plotCloud3D('reg','old',cols(i,:));
end
axis equal
title('Registered Coordinates')

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
for i = 1:length(newClouds)
    allVariances = [allVariances; newClouds(i).varRegXYZ];
    disp(['Mean Standard Deviation of X, Y, Z post registration (mm) for cloud number ', num2str(i), ':'])
    disp(mean(sqrt(newClouds(i).varRegXYZ))*1000)
end

disp('Mean Standard Deviation of X, Y, Z post registration for all clouds(mm):')
disp(mean(sqrt(allVariances))*1000)