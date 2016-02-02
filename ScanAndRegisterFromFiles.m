function newClouds = ScanAndRegisterFromFiles(varargin)
% newClouds = ScanAndRegisterFromFiles(scanfilename, planesfilename, spheresfilename [, homescanindex])
% This function takes in the file names of the scanners, planes and spheres
% to be used in a virtual scan and registration. It performs the scan,
% registers them together, and returns an array of Cloud objects which
% contain all of the spatial and statistical information of the clouds.
% Specifying a fourth input argument allows the user to specify the home
% scan index. If none is provided, then it is assumed to be 1. 

scanFileName = varargin{1};
planesFileName  = varargin{2};
spheresFileName = varargin{3};

if length(varargin) > 3
    homeScanIndex  = varargin{4};
else
    homeScanIndex = 1;
end

scannersRaw = load(scanFileName);

for i = 1:size(scannersRaw,1)
    scan(i) = Scanner(scannersRaw(i,1:3),scannersRaw(i,4:6));
    scan(i).angularIncrement = scan(i).angularIncrement*2;
    
    disp(['Scan number ',num2str(i)])
    tic
    clouds(i) = scan(i).ScanScene(planesFileName, spheresFileName);
    toc
end

[params, C_x] = registerScans(scan(homeScanIndex ), clouds, spheresFileName);

disp('Loading Data...')
scan = scan(1).updateScans(scan, params, C_x, homeScanIndex);

for i = 1:length(scan)
    if i ~= homeScanIndex
        newClouds(i) = Cloud(clouds(i).XYZ,scan(i),'xyz');
        newClouds(i).varXYZ = clouds(i).varXYZ;
        newClouds(i).varRTA = clouds(i).varRTA;
        newClouds(i).errXYZ = clouds(i).errXYZ;
        newClouds(i).errRTA = clouds(i).errRTA;
        [newClouds(i).varRegXYZ, newClouds(i).maxRegSTD] = newClouds(i).propRegErrors();
    else
        newClouds(i) = clouds(i);
        newClouds(i).regXYZ = newClouds(i).XYZ;
        newClouds(i).varRegXYZ = newClouds(i).varXYZ;
        newClouds(i).maxRegSTD = max(newClouds(i).varXYZ,[],2);
    end
    
end
disp('Process complete!')


end