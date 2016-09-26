classdef Cloud
    properties
        GLOBALXYZ;
        XYZ;
        regXYZ;
        RTA;
        varXYZ;
        varRTA;
        varRegXYZ;
        errXYZ;
        errRTA;
        scan;
        maxRegSTD;
    end
    
    methods
        function obj = Cloud(varargin)
            % This function creates a cloud object, which stores all the
            % imformation neccessary for the analysis of the performance in
            % a laser scanner. It stores the absolute XYZ positions of the
            % points (GLOBALXYZ), the xyz coordinatse of the points from
            % the perspective of the scanner (xyz), the range horz angle,
            % vert angle of the same points (RTA). It also stores
            % the variances of all those points. It stores the absolute
            % errors of XYZ and RTA as well, since those are only added as
            % a random number scaled to the expected variance in the rta
            % dimensions. This class also stores a scan, which contains all
            % the information about the scanner which created this cloud.
            % This, however, does not need to be filled, it is possible to
            % have a cloud with no information in the scan property. The
            % following input types are supported. The following
            % initialization are not case sensitive. All input types are to
            % be of a format of nx3, each row corresponding to a single
            % observation.
            %
            % cloud1 = Cloud(xyz, 'xyz'); cloud1 = Cloud(xyz, 'localxyz')
            % No global coordinates are calculated, and the scan
            % property is left blank, but the spherical coordinates are
            % calucaled.
            %
            % cloud1 = Cloud(rta, 'spherical'); cloud1 = Cloud(rta, 'angle')
            % Both of these two types behave similarly to the previous
            % ones, but the inputs are range, horizontal angel, vertical
            % angle. Again, no global coordinatse are calculated and the
            % scan property is left blank.
            %
            % cloud1 = Cloud(xyz, scan1, 'xyz')
            % This method now updates the scan properity with scan1, and
            % also calculates the global xyz coordinates based on that
            % scan's location and direction. This is supported for both the
            % xyz and the spherical coordinates as inputs.
            %
            % cloud1 = Cloud(globalxyz, scan1, 'globalxyz')
            % global xyz coordinates are valid as inputs only if a valid
            % scan is also included as input.
            %
            % If the scan provided also has registartion parameters, then
            % the coordinates of xyz after registartion are also
            % calculated.
            
            % Checking number of inputs
            if length(varargin) == 2
                input = varargin{1};
                inputType = varargin{2};
            elseif length(varargin) == 3
                input = varargin{1};
                TempScan= varargin{2};
                inputType = varargin{3};
            else
                error('Number of Input Arguments not supported')
            end
            
            % Checking the types of the input data
            if ~ischar(inputType)
                error('The second parameter, Input Type, must be a string');
            end
            
            if ~isnumeric(input)
                error('The first input must be numerical')
            end
            
            if size(input,2)~=3
                error('The input maxtrix must be size nx3')
            end
            
            if length(varargin) == 3
                if (isa(TempScan, 'Scanner'))
                    obj.scan = TempScan;
                else
                    error('If there are three elements, the second must be of type Scanner');
                    
                end
            end
            
            % Forcing the input string to be lowercase, so that its not
            % case sensitive
            inputType = lower(inputType);
            
            if strcmp(inputType,'angle') || strcmp(inputType,'spherical')
                % RTA is equal to the input
                obj.RTA = input;
                obj.XYZ = zeros(size(input));
                
                % XYZ is calculated from the input
                obj.XYZ(:,1) = input(:,1).*cos(input(:,2)).*cos(input(:,3));
                obj.XYZ(:,2) = input(:,1).*sin(input(:,2)).*cos(input(:,3));
                obj.XYZ(:,3) = input(:,1).*sin(input(:,3));
                
                % If a scan was provided, the global XYZ are also
                % calculated
                if length(varargin) == 3
                    
                    obj.GLOBALXYZ = zeros(size(input));
                    if size(input,1) ~= 0 & size(input,2)~=0
                        
                        tempXYZ = (rotz(-obj.scan.direction(3))*roty(-obj.scan.direction(2))...
                            *rotx(-obj.scan.direction(1))*(obj.XYZ'))';
                        
                        obj.GLOBALXYZ = [tempXYZ(:,1)+obj.scan.location(1),...
                            tempXYZ(:,2)+obj.scan.location(2),tempXYZ(:,3)+obj.scan.location(3)];
                    end
                end
                
            elseif strcmp(inputType,'xyz') || strcmp(inputType,'localxyz')
                % Input assigned to XYZ
                obj.XYZ = input;
                obj.RTA = zeros(size(input));
                
                % RTA calculated from RTA
                obj.RTA(:,1) = sqrt( input(:,1).^2 + input(:,2).^2 + input(:,3).^2 );
                obj.RTA(:,2) = atan2( input(:,2), input(:,1));
                obj.RTA(:,3) = atan2( input(:,3),  sqrt( input(:,1).^2 + input(:,2).^2));
                
                % If a scan was provided, calculate the global XYZ
                if length(varargin) == 3
                    
                    obj.GLOBALXYZ = zeros(size(input));
                    if size(input,1) ~= 0 & size(input, 2) ~= 0
                        tempXYZ = (rotz(-obj.scan.direction(3))*roty(-obj.scan.direction(2))...
                            *rotx(-obj.scan.direction(1))*(obj.XYZ'))';
                        
                        obj.GLOBALXYZ = [tempXYZ(:,1)+obj.scan.location(1),...
                            tempXYZ(:,2)+obj.scan.location(2),tempXYZ(:,3)+obj.scan.location(3)];
                    end
                end
                
            elseif strcmp(inputType,'globalxyz')
                obj.GLOBALXYZ = input;
                
                
                tempXYZ = [input(:,1) - obj.scan.location(1), input(:,2) - obj.scan.location(2),input(:,3) - obj.scan.location(3)];
                
                obj.XYZ = (rotz(obj.scan.direction(3))*roty(obj.scan.direction(2))...
                    *rotx(obj.scan.direction(1))*tempXYZ')';
                
                obj.RTA(:,1) = sqrt( obj.XYZ(:,1).^2 + obj.XYZ(:,2).^2 + obj.XYZ(:,3).^2 );
                obj.RTA(:,2) = atan2( obj.XYZ(:,2), obj.XYZ(:,1));
                obj.RTA(:,3) = atan2( obj.XYZ(:,3),  sqrt( obj.XYZ(:,1).^2 + obj.XYZ(:,2).^2));
                
                
            else
                error('Input type in not recognized, valid types are spherical, localxyz, or globalxyz');
            end
            
            % This may need some adjusting
            % If the provided scan also has registration parameters,
            % calculate the coordinates of the registered xyz points.
            if exist(obj.scan,'var') && obj.scan.regParams ~= zeros(6,1)
                obj.regXYZ = (rotz(obj.scan.regParams(6))*...
                    roty(obj.scan.regParams(5))*rotx(obj.scan.regParams(4))*...
                    (obj.XYZ') + repmat(obj.scan.regParams(1:3),1,size(obj.XYZ,1)))';
            end
            
        end
        
        function plotCloud3D(obj, varargin)
            % This is a fairly basic 3D plotting function. It simply plots
            % the x, y, and z coordinates in a new figure, and sets the
            % axis to equal. It accepts an input indicates which coordinate
            % system to use, local, global or registered. The second
            % possible input is to specify if the figure should be in a new
            % figure or not. Enter 'new', 'newfigure', or 'figure' for a
            % new figure, or 'old' to keep the old figure. Finally, the
            % colour can also be specified in the final input argument. It
            % should be a 1x3 row vector scaled between 0 and 1.
            
            if size(varargin,1)>0
                if ischar(varargin{1})
                    input = varargin{1};
                    input = lower(input);
                else
                    warning('Input type not recognized, defaulting to local coordinates');
                    input = 'xyz';
                end
                
                if strcmp(input ,'xyz') | strcmp(input,'local') | strcmp(input,'localxyz')
                    data = obj.XYZ;
                    scanLoc = [0,0,0];
                elseif strcmp(input ,'globalxyz') | strcmp(input,'global')
                    data = obj.GLOBALXYZ;
                    scanLoc = obj.scan.location;
                elseif strcmp(input, 'reg') | strcmp(input, 'regxyz') |strcmp(input, 'registered')
                    data = obj.regXYZ;
                    scanLoc = obj.scan.regParams(1:3);
                else
                    warning('Input format not recognized, defaulting to local coordinates');
                    data = obj.XYZ;
                    scanLoc = [0,0,0];
                end
            else
                input = 'xyz';
                data = obj.XYZ;
            end
            
            if length(varargin)>1
                if strcmp(varargin{2},'new')|strcmp(varargin{2},'newfigure')|strcmp(varargin{2},'figure')
                    figure,
                end
            end
            
            
            if length(varargin)>2
                if size(varargin) == [1,3];
                    col = varargin{3};
                end
            else
                col = [0,0,1];
            end
            
            plot3(data(:,1),data(:,2), data(:,3), '.','Color',col)
            plot3(scanLoc(1),scanLoc(2),scanLoc(3),'<','Color',col)
            
            
            if length(varargin)>1
                if strcmp(varargin{2},'new')|strcmp(varargin{2},'newfigure')|strcmp(varargin{2},'figure')
                    axis equal
                    xlabel('X axis (m)')
                    ylabel('Y axis (m)')
                    zlabel('Z axis (m)')
                    title(input)
                end
            end
        end
        
        function plotCloud2D(obj, varargin)
            % This function creates a 2D plot of the point cloud, with
            % horizontal angle along the x axis and vertical angle angle
            % the y axis. The each point gets a colour, which corresponds
            % (in a non-linear way) to the range of the point from the
            % scanner.
            
            % Here you can specify the number of colours used in the plot.
            % A lower number will have less detail, while a higher number
            % will take slightly more computation.
            if size(varargin,1) >= 1
                numColours = varargin{1};
            else
                numColours = 100;
            end
            
            temp = sortrows(obj.RTA(:,1));
            
            colInc = round(length(temp)/numColours);
            
            colours = zeros(numColours, 1);
            for i = 1:numColours
                if i*colInc <= size(temp,1)
                    colours(i) = temp(colInc*i);
                else
                    colours(i) = temp(end);
                    break
                end
            end
            
            figure,
            hold on
            for i=1:numColours
                if i == 1
                    key =  obj.RTA(:,1) < colours(i);
                elseif i == numColours
                    key = obj.RTA(:,1) >= colours(i-1);
                else
                    key = (obj.RTA(:,1) < colours(i)) .* (obj.RTA(:,1) >= colours(i-1));
                end
                
                data = deleteRowKey([obj.RTA(:,2), obj.RTA(:,3)] , key);
                
                plot(data(:,1),  data(:,2), '.', 'Color', [sqrt(1 - ((i-1)/(numColours-1))), 0,  0]);
            end
            xlabel('Horizontal Angle (radians)')
            ylabel('Vertical Angle (radians)')
            %             set(gca,'Xlim',[-pi,pi]);
            %             set(gca,'YLim',obj.scan.verticalRange);
            
        end
        
        function newCloud = combineClouds(obj, cloud2)
            % The combine clouds method takes a cloud as input. That cloud
            % must have a scan value equal to that of the cloud calling
            % this function. This function then returns a new cloud which
            % combines the two clouds, having removed any point which has
            % the same horizontal and vertical value as a point already
            % present. The point that is removed is the point with the
            % larger range. Unfortunately, it is not possible to simply add
            % the second cloud to the first. They must be returned by this
            % function.
            
            if sum((obj.scan.location == cloud2.scan.location).*(obj.scan.direction == cloud2.scan.direction)) ~= 3
                error('Scans must be equivalent');
            end
            
            % Combining Clouds into vectors
            if size(obj.errRTA,1) ~= 0
                cloudsAngular = [obj.RTA-obj.errRTA; cloud2.RTA - cloud2.errRTA];
            else
                cloudsAngular = [obj.RTA; cloud2.RTA];
            end
            
            varRTA = [obj.varRTA; cloud2.varRTA];
            varXYZ = [obj.varXYZ; cloud2.varXYZ];
            
            errRTA = [obj.errRTA; cloud2.errRTA];
            errXYZ = [obj.errXYZ; cloud2.errXYZ];
            
            % Sorting All clouds according to horizontal angle, vertical
            % angle, then range.
            DiscreteCloudsAngular = round([cloudsAngular(:,1)*10000,...
                cloudsAngular(:,2)/obj.scan.angularIncrement,  cloudsAngular(:,3)/obj.scan.angularIncrement]);
            
            allSortedData = real(sortrows([DiscreteCloudsAngular,varRTA,varXYZ,errRTA,errXYZ],[2,3,1]));
            
            % Separating data
            sortedCloudsAngular = allSortedData (:,1:3);
            sortedVarRTA = allSortedData (:,4:6);
            sortedVarXYZ = allSortedData (:,7:9);
            sortedErrRTA = allSortedData (:,10:12);
            sortedErrXYZ = allSortedData (:,13:15);
            
            % Identifying elements with matching horizontal and vertical
            % angles
            matchingAngles1 = (sortedCloudsAngular(1:end-1,3) ~= sortedCloudsAngular(2:end,3))|...
                (sortedCloudsAngular(1:end-1,2) ~= sortedCloudsAngular(2:end,2));
            
            % Deleting duplicate elements
            matchingAngles = [1; matchingAngles1];
            newDiscreteCloudAng = deleteRowKey(sortedCloudsAngular, matchingAngles);
            newVarRTA = deleteRowKey(sortedVarRTA, matchingAngles);
            newVarXYZ = deleteRowKey(sortedVarXYZ, matchingAngles);
            newErrRTA = deleteRowKey(sortedErrRTA, matchingAngles);
            newErrXYZ = deleteRowKey(sortedErrXYZ, matchingAngles);
            
            newCloudAng = [newDiscreteCloudAng(:,1)/10000,...
                newDiscreteCloudAng(:,2)*obj.scan.angularIncrement,...
                newDiscreteCloudAng(:,3)*obj.scan.angularIncrement];
            newCloudAng = newCloudAng+newErrRTA;
            
            % Storing non-duplicate elements in new cloud
            newCloud = Cloud(newCloudAng,obj.scan,'angle');
            newCloud.varRTA = newVarRTA;
            newCloud.varXYZ = newVarXYZ;
            newCloud.errRTA = newErrRTA;
            newCloud.errXYZ = newErrXYZ;
            
            if ~isreal(newErrXYZ)
                print('weird')
            end
            
        end
        
        function [noiseRTA, noiseXYZ] = createNoise(obj)
            % create noise is a function that takes the variance of the
            % range, horizontal and vertical angles and returns a random
            % sample of noise which follows the variances of those
            % properties, then converts those errors into xyz coordinates.
            noiseRTA = randn(size(obj.varRTA)).*sqrt(abs(obj.varRTA));
            RTA = obj.RTA;
            
            noiseXYZ = zeros(size(noiseRTA));
            
            a =  zeros(3);
            
            for i = 1:size(noiseRTA,1)
                
                L = noiseRTA(i,:)';
                
                a(1,1) = cos(RTA(i,3))*cos(RTA(i,2));
                a(2,1) = cos(RTA(i,3))*sin(RTA(i,2));
                a(3,1) = sin(RTA(i,3));
                
                a(1,2) = -RTA(i,1)* cos(RTA(i,3))*sin(RTA(i,2));
                a(2,2) = RTA(i,1)* cos(RTA(i,3))*cos(RTA(i,2));
                
                a(1,3) = -RTA(i,1)*sin(RTA(i,3))*cos(RTA(i,3));
                a(2,3) = -RTA(i,1)*sin(RTA(i,3))*sin(RTA(i,3));
                a(3,3) = RTA(i,1)*cos(RTA(i,3));
                
                noiseXYZ(i,:) = (a*L)';
                
            end
            
        end
        
        function newCloud = deletePoints(obj,key)
            % Delete points is a method that deletes a point from a cloud,
            % based on a logical index the length of the cloud. If an index
            % is 1, the point in kept. If it is zero, the point is deleted.
            %
            if size(obj.XYZ,1) ~= size(key,1)
                error('Key must be the same length as the Cloud');
            end
            
            allData = [obj.XYZ,obj.varRTA,obj.varXYZ,obj.errRTA,obj.errXYZ];
            
            newData = deleteRowKey(allData, key);
            
            % Storing non-duplicate elements in new cloud
            newCloud = Cloud(newData(:,1:3),obj.scan,'xyz');
            newCloud.varRTA = newData(:,4:6);
            newCloud.varXYZ = newData(:,7:9);
            newCloud.errRTA = newData(:,10:12);
            newCloud.errXYZ = newData(:,13:15);
            
        end
        
        function sphereCloud = collectSpherePoints(obj,sphere)
            % Given a sphere object, this function collects all of the
            % point that fall within a range from the center of the sphere
            % up to the radius, plus 10% to accomodate for noise.
            sphereCenterCloud = Cloud(sphere.center,obj.scan,'GlobalXYZ');
            xyz = sphereCenterCloud.XYZ;
            
            dist = sqrt( (obj.XYZ(:,1)-xyz(1)).^2 + (obj.XYZ(:,2)-xyz(2)).^2 + (obj.XYZ(:,3)-xyz(3)).^2 );
            
            key = (dist/1.1)<sphere.radius;
            
            sphereCloud = obj.deletePoints(key);
        end
        
        function [center, variance] = modelSphere(obj,sphere)
            % Assuming the cloud is spherical in nature, given an input
            % sphere class, the center of the cloud according to the radius
            % of the sphere class is estimated, along with the variance of
            % the xyz position.
            
            xyz = obj.XYZ;
            sphereCenterCloud = Cloud(sphere.center,obj.scan,'GlobalXYZ');
            TrueCenter = sphereCenterCloud.XYZ;
            
            % contrain the radius to the known value
            rcon = sphere.radius;
            
            % 1mm std for the radius of targets
            sd = 0.001;
            
            n_pt_obs=max(size(xyz));
            
            if (n_pt_obs <= 3)
                converged=0;
                x0=nan(4,1);
                center = x0(1:3);
                variance = center;
                disp('Sphere Fit: Insufficient number of Points')
                return
            end
            
            % circle fit
            tol=1e-4;
            max_iter=12;
            
            xc= TrueCenter(1) + (randn*sphere.radius*0.05);
            yc= TrueCenter(2)+ (randn*sphere.radius*0.05);
            zc= TrueCenter(3)+ (randn*sphere.radius*0.05);
            
            if rcon > 0
                r=rcon;
                wt=1/sd^2;
            else
                r=max(range(xyz))/2;
            end
            
            % order: xc yc zc r
            C_L = zeros(size(obj.XYZ,1));
            x0=[ xc yc zc r ]';
            
            % This calculation of the observation variance may need some
            % theoretical backing up.
            
            %             for i = 1:size(obj.XYZ,1);
            %                 C_L(i,i) = sum(obj.varXYZ(i,:));
            %             end
            
            C_L = diag(sum(obj.varXYZ,2));
            C_L_inv = pinv(C_L);
            
            converged=0;
            for i=1:max_iter
                A=zeros(n_pt_obs,4);
                w=zeros(n_pt_obs,1);
                
                xc=x0(1);
                yc=x0(2);
                zc=x0(3);
                r=x0(4);
                
                for j=1:n_pt_obs
                    x=xyz(j,1);
                    y=xyz(j,2);
                    z=xyz(j,3);
                    
                    dx=x-xc;
                    dy=y-yc;
                    dz=z-zc;
                    
                    A(j,1)=-2*dx;
                    A(j,2)=-2*dy;
                    A(j,3)=-2*dz;
                    A(j,4)=-2*r;
                    
                    w(j)=dx^2+dy^2+dz^2-r^2;
                end
                
                N=A'*C_L_inv*A;
                U=A'*C_L_inv*w;
                % add radius constraint
                if (rcon > 0)
                    N(4,4)=N(4,4)+wt;
                    U(4,1)=U(4,1)+wt*(r-rcon);
                end
                
                % check for singularity
                rc=rcond(N);
                if ( (rc < eps) || isnan(rc) )
                    return
                end
                %[i  cond(A'*A) rcond(A'*A)]
                
                dx=-inv(N)*U;
                x0=x0+dx;
                
                %[i dx']
                
                if (max(abs(dx)) < tol)
                    converged=1;
                    break
                end
            end
            center = x0(1:3);
            variance = diag(pinv(N(1:3,1:3)));
            
            if converged == 0
                error('Sphere Fit: Failed to converge')
            end
        end
        
        function [center, variance] = modelSphereAdvanced(obj,sphere)
            % Assuming the cloud is spherical in nature, given an input
            % sphere class, the center of the cloud according to the radius
            % of the sphere class is estimated, along with the variance of
            % the xyz position.
            % This 'Advanced' version of the sphere model function attempts
            % to be more rigorous by including the co-variances between
            % the object point coordinates. 
            
            varCovar = VarianceCalculationSphereCell(obj.scan, obj, sphere);
            xyz = obj.XYZ;
            
            sphereCenterCloud = Cloud(sphere.center,obj.scan,'GlobalXYZ');
            TrueCenter = sphereCenterCloud.XYZ;
            
            
            % contrain the radius to the known value
            rcon = sphere.radius;
            
            % 1mm std for the radius of targets
            sd = 0.001;
            
            n_pt_obs=max(size(xyz));
            
            if (n_pt_obs <= 3)
                converged=0;
                x0=nan(4,1);
                center = x0(1:3);
                variance = center;
                disp('Sphere Fit: Insufficient number of Points')
                return
            end
            
            % circle fit
            tol=1e-4;
            max_iter=12;
            
            xc= TrueCenter(1) + (randn*sphere.radius*0.05);
            yc= TrueCenter(2)+ (randn*sphere.radius*0.05);
            zc= TrueCenter(3)+ (randn*sphere.radius*0.05);
            
            if rcon > 0
                r=rcon;
                wt=1/sd^2;
            else
                r=max(range(xyz))/2;
            end
            
            % order: xc yc zc r
            x0=[ xc yc zc r ]';
            
            % Using varCovar
            C_L = [];
            for i = 1:length(varCovar)
                C_L = blkdiag(C_L, varCovar{i});
            end
            
            C_L_inv = pinv(C_L);
            
            converged=0;
            for i=1:max_iter
%                 A=zeros(n_pt_obs*3,4);
%                 w=zeros(n_pt_obs*3,1);
                A = [];
                w = [];
                
                xc=x0(1);
                yc=x0(2);
                zc=x0(3);
                r=x0(4);
                
                for j=1:n_pt_obs
                    subA = zeros(3,4);
                    subw = zeros(3,1);
                    
                    x=xyz(j,1);
                    y=xyz(j,2);
                    z=xyz(j,3);
                    
                    dx=x-xc;
                    dy=y-yc;
                    dz=z-zc;
                    
                    dnom = sqrt(r^2 - dy^2 - dz^2);
                    subA(1,1) = 1;
                    subA(1,2) = dy/dnom;
                    subA(1,3) = dz/dnom;
                    subA(1,4) =  r/dnom;
                    
                    subw(1) = dnom+xc;
                    
                    dnom = sqrt(r^2 - dx^2 - dz^2);
                    subA(2,1) = dx/dnom;
                    subA(2,2) = 1;
                    subA(2,3) = dz/dnom;
                    subA(2,4) =  r/dnom;
                    
                    subw(2) = dnom+yc;
                    
                    dnom = sqrt(r^2 - dx^2 - dy^2);
                    subA(3,1) = dx/dnom;
                    subA(3,2) = dy/dnom;
                    subA(3,3) = 1;
                    subA(3,4) =  r/dnom;
                    
                    subw(3) = dnom+zc;
                    
                    A = [A;subA]
                    
                    w = [w;subw];
                end
                
                N=A'*C_L_inv*A;
                U=A'*C_L_inv*w;
                % add radius constraint
                if (rcon > 0)
                    N(4,4)=N(4,4)+wt;
                    U(4,1)=U(4,1)+wt*(r-rcon);
                end
                
                % check for singularity
                rc=rcond(N);
                if ( (rc < eps) || isnan(rc) )
                    return
                end
                %[i  cond(A'*A) rcond(A'*A)]
                
                dx=-inv(N)*U;
                x0=x0+dx;
                
                %[i dx']
                
                if (max(abs(dx)) < tol)
                    converged=1;
                    break
                end
            end
            center = x0(1:3);
            variance = diag(pinv(N(1:3,1:3)));
            
            if converged == 0
                error('Sphere Fit: Failed to converge')
            end
        end
        
        function [estSphereCenters, sphereVariances] = estimateSphereTargets(obj, AllSpheres)
            % This function takes an array of sphere objects, and esitmates
            % the center and associated variances for each one.
            numSpheres = length(AllSpheres);
            estSphereCenters = zeros(numSpheres,3);
            sphereVariances = estSphereCenters;
            
            for i = 1:numSpheres
                sphereCloud = obj.collectSpherePoints(AllSpheres(i));
                [estSphereCenters(i,:),sphereVariances(i,:)] = sphereCloud.modelSphere(AllSpheres(i));
            end
            
        end
        
        function [varRegXYZ, maxRegSTD] = propRegErrors(obj)
            % This function takes the xyz coordinates of a cloud, as well
            % as the registration parameters and associated variances and
            % calculates the variances of the xyz coordinates after
            % registration.
            varRegXYZ = zeros(size(obj.XYZ));
            maxRegSTD = zeros(size(obj.XYZ,1),1);
            
            for i = 1:length(obj.XYZ)
                
                X = obj.XYZ(i,1);
                Y = obj.XYZ(i,2);
                Z = obj.XYZ(i,3);
                
                w = obj.scan.regParams(4);
                ph = obj.scan.regParams(5);
                k = obj.scan.regParams(6);
                
                A = zeros(3,9);
                
                A(1:3, 1:3) = eye(3);
                
                A(1,4) = (Z*(cos(w)*sin(k) + cos(k)*sin(w)*sin(ph))) - (Y*(sin(k)*sin(w) - cos(k)*cos(w)*sin(ph)));
                A(1,5) = (Y*(cos(k)*cos(ph)*sin(w))) - (Z*(cos(k)*cos(w)*cos(ph)))- (X*(cos(k)*sin(ph)));
                A(1,6) = (Y*(cos(k)*cos(w)-sin(k)*sin(w)*sin(ph))) + (Z*(cos(k)*sin(w)+cos(w)*sin(k)*sin(ph))) - (X*(cos(ph)*sin(k)));
                
                A(2,4) = (Z*(cos(k)*cos(w) -sin(k)*sin(w)*sin(ph))) -(Y*(cos(k)*sin(w) + cos(w)*sin(k)*sin(ph)));
                A(2,5) = (X*(sin(k)*sin(ph))) + (Z*(cos(w)*cos(ph)*sin(k))) + (Y*(cos(ph)*sin(k)*sin(w)));
                A(2,6) = -(Y*(cos(w)*sin(k) + cos(k)*sin(w)*sin(ph))) - (Z*(sin(k)*sin(w)-cos(k)*cos(w)*sin(ph))) - (X*(cos(k)*cos(ph)));
                
                A(3,4) = -(Y*(cos(w)*cos(ph))) - (Z*(cos(ph)*sin(w)));
                A(3,5) = (X*cos(ph)) - (Z*(cos(w)*sin(ph))) + (Y*(sin(w)*sin(ph)));
                
                A(1:3,7:9) = rotz(k)*roty(ph)*rotx(w);
                
                C_L = blkdiag(obj.scan.regParamVarCovar,diag( obj.varXYZ(i,:)));
                
                if rcond(C_L) < 1e-12
                    disp(['Poor condition on C_L matrix for scan point number ', num2str(i)])
                end
                
                C_x = A * C_L * A';
                
                varRegXYZ(i,:) = diag(C_x)';
                
                maxVar3 = eig(C_x);
                maxRegSTD(i) = sqrt(maxVar3(1));
                
            end
            
        end
    end
    
    methods (Static)
        function varXYZ = getXYZvariance(obj)
            % This function calculates the xyz variances based on the
            % variances of range, horizontal and vertical angle.
            RTA = obj.RTA;
            varRTA = obj.varRTA;
            
            errors = zeros(size(varRTA));
            
            a =  zeros(3);
            
            for i = 1:size(varRTA,1)
                
                C_L = diag(varRTA(i,:));
                
                a(1,1) = cos(RTA(i,3))*cos(RTA(i,2));
                a(2,1) = cos(RTA(i,3))*sin(RTA(i,2));
                a(3,1) = sin(RTA(i,3));
                
                a(1,2) = -RTA(i,1)* cos(RTA(i,3))*sin(RTA(i,2));
                a(2,2) = RTA(i,1)* cos(RTA(i,3))*cos(RTA(i,2));
                
                a(1,3) = -RTA(i,1)*sin(RTA(i,3))*cos(RTA(i,3));
                a(2,3) = -RTA(i,1)*sin(RTA(i,3))*sin(RTA(i,3));
                a(3,3) = RTA(i,1)*cos(RTA(i,3));
                
                C_x = (a'* C_L^-1 *a)^-1;
                
                errors(i,:) = diag(C_x)';
            end
            
            varXYZ = errors;
        end
        
        
        
    end
end