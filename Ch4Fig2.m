% chapter 4 figure 2: TLS observation coordinate system
sphere = Sphere([2,0,0],0.076)

scanner = Scanner([0,0,0],[0,0,0])

cloud = scanner.ScanSphere(sphere)

figure, subplot(1,1,1,'FontSize',20)
plot3(cloud.XYZ(:,1),cloud.XYZ(:,2),cloud.XYZ(:,3), '.')
axis equal