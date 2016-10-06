% Display Scene

close all
clear
clc


%% load the scene
scansraw = load('kuukpak scan models.txt');
planesraw = load('kuukpak planar model.txt');
spheresraw = load('kuukpak sphere models.txt');

%% display scene elements

figure, subplot(1,1,1,'FontSize',16), hold on

% Display the planes (crude, but it works)
plot(planesraw(:,1), planesraw(:,2), 'LineWidth',2)

% Display the Scanner Locations
plot(scansraw(:,1), scansraw(:,2), 'r<','LineWidth',2)

% Display Spheres
plot(spheresraw(:,1), spheresraw(:,2), 'ko','LineWidth',2)

% Format images nicely 
axis equal
legend('Plane Boundaries','Scan Locations','Spherical Targets','Location','EastOutside')
xlabel('X (m)')
ylabel('Y (m)')
title('Layout of Simulated Scene')
