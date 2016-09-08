classdef Scanner
    properties
        verticalRange
        horizontalRange
        maxDistance
        angularIncrement
        location
        direction
        obsStd
        levelStd
        regParams
        regParamVarCovar
        
    end
    methods
        function obj = Scanner(loc,dir)
            % The Scanner class contains all relevant information about the
            % scan. This includes the location (xyz) of the perpective
            % center in global coordinates, the direction (omega, phi,
            % kappa) of the scanning reference axis. The system assumes
            % that the scanner is level, so omega and phi are always zero.
            % Also stored and the vertical and horizontal ranges, the
            % maximum range, the angular increment, the standard deviations
            % of the range, horz and vert angles, the standard devitation
            % of the leveling device, and well as the registration
            % parameters and associated variances.
            if nargin > 1
                if isnumeric(loc) &&  isnumeric(dir)
                    obj.location = loc;
                    obj.direction = dir*pi/180; % Angle are entered in degrees, and converted to rads
                    obj.verticalRange = [-27.5,90]*pi/180;
                    obj.horizontalRange = [-pi,pi];
                    obj.maxDistance = 120;
                    obj.angularIncrement = 0.009*pi/180 * 4; % Less thas this and things start to crash
                    obj.obsStd = [1.7e-3, 3.2e-4, 2.2e-4]; % meters, radians, radians
                    obj.levelStd = [1.8e-4, 1.8e-4, 0];
                    obj.regParams = zeros(6,1);
                    obj.regParamVarCovar = zeros(6);
                    
                else
                    error('Location and Direction must be numeric')
                end
            end
        end
        
        function cloud = ScanPlane(obj, plane)
            % Scan plane takes a plane object, and produces a cloud object
            % of it. This version is an attempt at being more robust
            
            % Check normal direction
            if dot( plane.centroid - obj.location, plane.normal)>0
                cloud = Cloud(zeros(0,3),obj,'angle');
                return
            end
            
            % Check plane distance i.e. closest point on plane < max range
            n = plane.normal;
            S = obj.location;
            min_distance_to_scanner = abs(abs(dot(n,S)) - plane.distance_to_origin);
            
            if min_distance_to_scanner > obj.maxDistance
                cloud = Cloud(zeros(0,3),obj,'angle');
                return
            end
            
            % Delete angular increments outside the extremes (may be too
            % much work, since edges may be out, while vertexes may be in)
            
            % Subtract the rotation of the scanner to account for change in
            % direction between the scanner and the object space coords
            
            % Also, these minimum angular values need to be rounded to the
            % nearest increment
            
            % This is the points of interest method. The lines between the
            % vertexes of the plane are checked for extremes in horizontal
            % and vertical angle. These determine which points we check to
            % see if the fall within the plane.
            
            newVertexs = Cloud(plane.vertexs, obj, 'GlobalXYZ');
            newPlane = Plane(newVertexs.XYZ);
            n = newPlane.normal;
            
            V = newVertexs.XYZ;
            
            v1 = repmat(V(1,:),1000,1);
            v2 = repmat(V(2,:),1000,1);
            
            v12 = repmat(V(2,:) - V(1,:),1000,1);
            t = [0.001:0.001:1]';
            t3 = [t,t,t];
            line12 = v1+(v12.*t3);
            
            v13 = repmat(V(3,:) - V(1,:),1000,1);
            line13 = v1+(v13.*t3);
            
            v23 = repmat(V(3,:) - V(2,:),1000,1);
            line23 = v2+(v23.*t3);
            
            poi = [V;line12;line13;line23];
            
            hAngle = atan2(poi(:,2), poi(:,1));
            h_min = min(hAngle);
            h_min = ceil(h_min/obj.angularIncrement)*obj.angularIncrement;
            
            h_max = max(hAngle);
            h_max = floor(h_max/obj.angularIncrement)*obj.angularIncrement;
            
            if (newPlane.centroid(1))<0
                hAngle = ((hAngle < 0)*2*pi)+hAngle;
                h_min = min(hAngle);
                h_max = max(hAngle);
            end
            
            vAngle = atan2(poi(:,3), sqrt((poi(:,1)).^2 + (poi(:,2)).^2 ));
            v_min = min(vAngle);
            % This is to ensure the plane is not scanned below the minimum
            % angle of the scanner.
            v_min = max([v_min, obj.verticalRange(1)]);
            
            v_max = max(vAngle);
            
            v_min_sign = sign(v_min);
            v_max_sign = sign(v_max);
            
            % Rounding the minimum and maximum coordinates to nearest
            % appropriate angular incements
            v_min = v_min_sign*floor(abs(v_min)/obj.angularIncrement)*obj.angularIncrement;
            
            v_max = v_max_sign*(floor(abs(v_max)/obj.angularIncrement)*obj.angularIncrement);
            
            % Initialize point cloud, in angular units
            vert = [v_min:obj.angularIncrement:v_max];
            horz = [h_min:obj.angularIncrement:h_max];
            
            vertical = repmat(vert', 1, length(horz));
            horizontal = repmat(horz, length(vert), 1);
            range = ones(size(horizontal))*(obj.maxDistance + 1);
            
            cloudAng = [range(:), horizontal(:), vertical(:)];
            
            % Determine the range of a ray to the plane.
            numerator =  newPlane.distance_to_origin;
            
            R = [cos(cloudAng(:,2)).*cos(cloudAng(:,3)),...
                sin(cloudAng(:,2)).*cos(cloudAng(:,3)),...
                sin(cloudAng(:,3))];
            
            denom = dot(repmat(n,size(R,1),1),R,2);
            
            t = - repmat(numerator,size(denom,1),1)./denom;
            
            % Determining coordinates of the intersection of ray and plane.
            new_point = [t.*(cos(cloudAng(:,2)).*cos(cloudAng(:,3))),...
                t.*(sin(cloudAng(:,2)).*cos(cloudAng(:,3))),...
                t.*sin(cloudAng(:,3))];
            
            % Using Barycentric Technique to determine if internal or
            % external to the vertixes
            
            v0 = repmat(V(3,:) - V(1,:), size(new_point,1),1);
            v1 = repmat(V(2,:) - V(1,:), size(new_point,1),1);
            v2 = new_point - repmat(V(1,:), size(new_point,1),1);
            dot00 = dot(v0,v0,2);
            dot01 = dot(v0,v1,2);
            dot11 = dot(v1,v1,2);
            dot02 = dot(v0,v2,2);
            dot12 = dot(v1,v2,2);
            
            invDenom = 1./((dot00.*dot11) - (dot01.*dot01));
            u = ((dot11.*dot02)-(dot01.*dot12)).*invDenom;
            v = ((dot00.*dot12)-(dot01.*dot02)).*invDenom;
            
            for i = 1:size(cloudAng,1)
                if (u(i)>0 & v(i)>0 & (u(i)+v(i) < 1))
                    cloudAng(i,1) = t(i);
                end
            end
            
            cloudAngFinal = deleteRowKey(cloudAng, (cloudAng(:,1)<obj.maxDistance)&(cloudAng(:,1)>0));
            cloudAngFinal(:,2:3) = obj.wrapToPi(cloudAngFinal(:,2:3));
            
            cloud = Cloud(cloudAngFinal,obj,'angle');
            
            [cloud.varRTA, cloud.varXYZ] = VarianceCalculationPlane(obj, cloud, plane);
            
            [cloud.errRTA, cloud.errXYZ] = cloud.createNoise();
            
            cloud.RTA = cloud.RTA + cloud.errRTA;
            cloud.XYZ= cloud.XYZ+ cloud.errXYZ;
        end
        
        function [varRTA, varXYZ]  = VarianceCalculationPlane(obj, cloud, plane)
            % The variance of the points on a plane are calculated here.
            % This is neccessary because the angle of incidence impacts the
            % variance in the range direction. Since the angle of incidence
            % on the surface of a plane and on the surface of a sphere are
            % calculated differently, there are two different functions for
            % their calculation.
            
            RTA = cloud.RTA;
            
            % The normal angle of the plane to the scanner's coordinate
            % system.
            newVertexs = Cloud(plane.vertexs,obj,'GlobalXYZ');
            newPlane = Plane(newVertexs.XYZ);
            n = newPlane.normal;
            
            
            % Calculate angle of incidence
%             incidence = zeros(size(RTA,1),1);
            
            %             for i = 1:length(incidence)
            %                 incidence(i) = acos(dot(n, -[cos(RTA(i,2)).*cos(RTA(i,3)),sin(RTA(i,2)).*cos(RTA(i,3)), sin(RTA(i,3))]));
            %             end
            
            incidence = acos( dot(repmat(n,size(RTA,1),1),-[cos(RTA(:,2)).*cos(RTA(:,3)),sin(RTA(:,2)).*cos(RTA(:,3)), sin(RTA(:,3))],2));
            
            % Calculate variance of alpha and theta
            varTheta = obj.obsStd(2).^2*ones(size(incidence));
            varAlpha = obj.obsStd(3).^2*ones(size(incidence));
            
            % The variance in range is scaled by the secant of the angle of
            % incidence, squared
            varR = (obj.obsStd(1) .*sec(incidence)).^2;
            
            varRTA = [varR, varTheta, varAlpha];
            
            varXYZ = zeros(size(varRTA));
            
            % A is the design matrix.
            a =  zeros(3);
            
            % This part is slow, consider performance enhancing alternative
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
                
                C_x = a* C_L *a';
                
                varXYZ(i,:) = diag(C_x)';
            end
            
            
        end
        
        function cloud = ScanSphere(obj, sphere)
            % This method takes in a sphere object and returns a cloud of
            % coordinates that are incident on its surface.
            S = obj.location;
            
            % Check distance to center - radius < max range
            if norm(sphere.center - S) - sphere.radius > obj.maxDistance
                cloud = Cloud(zeros(0,3),obj,'xyz');
                return
            end
            
            % Determine angular footprint of the target
            SFP = atan2(sphere.radius,norm(sphere.center - S));
            
            temp = Cloud(sphere.center,obj,'globalxyz');
            
            SCang = temp.RTA;
            inc = obj.angularIncrement;
            
            % Determine nearest valid angular increment in the vertical
            % directions
            Vmin = ceil((SCang(3) - SFP)/inc)*inc;
            
            Vmax = floor((SCang(3) + SFP)/inc)*inc;
            
            V_inc = (Vmin:inc:Vmax) -  SCang(3);
            %             V_inc = (Vmin:inc:Vmax);
            
            % Determining nearest horizontal angular increment, considering
            % the point density changes with the vertical angle.
            diff_H_ang =  max(acos(cos(V_inc) * cos(SFP))./cos(V_inc + SCang(3)));
            Hmin= ceil((SCang(2) - diff_H_ang)/inc)*inc;
            Hmax = floor((SCang(2) + diff_H_ang)/inc)*inc;
            
            
            % Declare point cloud, in angular units
            
            Vmin = ceil( Vmin/obj.angularIncrement )*obj.angularIncrement;
            
            Vmax = floor( Vmax/obj.angularIncrement )*obj.angularIncrement;
            
            Hmin = ceil( Hmin/obj.angularIncrement )*obj.angularIncrement;
            Hmax = floor( Hmax/obj.angularIncrement )*obj.angularIncrement;
            
            vert = [Vmin:obj.angularIncrement:Vmax];
            horz = [Hmin:obj.angularIncrement:Hmax];
            
            vertical = repmat(vert', 1, length(horz));
            horizontal = repmat(horz, length(vert), 1);
            range = ones(size(horizontal))*(obj.maxDistance + 1);
            
            cloudAng = [range(:), horizontal(:), vertical(:)];
            
            H = cloudAng(:,2) - SCang(2);
            V = cloudAng(:,3) - SCang(3);
            %             H = cloudAng(:,2);
            %             V = cloudAng(:,3);
            
            % delete points that have a greater angle between them, the
            % center of the scanner and the center of the sphere.
            % i.e. All points on the surface of the sphere should have an
            % angle less than the footprint of the sphere.
            cloud_ang_dist = acos( cos(V) .* cos(H.*cos(cloudAng(:,3))));
            
            points_on_sphere = deleteRowKey(cloudAng, cloud_ang_dist < SFP);
            
            % for all remaining angular increments, determine the point
            % range using the sphere radius
            sphere_angular_dist = deleteRowKey( cloud_ang_dist, cloud_ang_dist < SFP);
            
            beta = asin(SCang(1) * sin(sphere_angular_dist)/sphere.radius);
            beta2  = pi - beta;
            
            interior_angle = pi-beta-sphere_angular_dist;
            interior_angle2 = pi-beta2-sphere_angular_dist;
            
            distance_to_sphere_surface = sphere.radius .* sin(interior_angle)...
                ./sin(sphere_angular_dist);
            
            distance_to_sphere_surface2 = sphere.radius .* sin(interior_angle2)...
                ./sin(sphere_angular_dist);
            
            distance_true = min( [distance_to_sphere_surface, distance_to_sphere_surface2], [], 2);
            
            for i = 1:length(sphere_angular_dist)
                if abs(sin(sphere_angular_dist(i)))<1e-15
                    distance_true(i) = norm(sphere.center - S) - sphere.radius;
                end
            end
            
            cloudAng = real([distance_true, points_on_sphere(:,2:3)]);
            cloudAng = deleteRowKey(cloudAng, cloudAng(:,1)<obj.maxDistance);
            
            cloudAng(:,2:3) = obj.wrapToPi(cloudAng(:,2:3));
            
            cloud = Cloud(cloudAng, obj, 'angle');
            
            [cloud.varRTA, cloud.varXYZ] = VarianceCalculationSphere(obj, cloud, sphere);
            
            [cloud.errRTA, cloud.errXYZ] = cloud.createNoise();
            
            cloud.RTA = cloud.RTA + cloud.errRTA;
            cloud.XYZ= cloud.XYZ+ cloud.errXYZ;
        end
        
        function [varRTA, varXYZ] = VarianceCalculationSphere(obj, cloud, sphere)
            % This process calculates the variance of the points on the
            % surface of a sphere. Since the angle of incidence of the
            % laser on the surface being scanned plays a part in
            % determining the variance in the range direction, and the
            % angle of incidence is calculated on different for a sphere
            % than a plane, two different functions are needed.
            
            % Converts the coordinates of the sphere into local
            % coordinates,
            SCcloud = Cloud(sphere.center,obj,'GlobalXYZ');
            SC = SCcloud.XYZ;
            
            RTA = cloud.RTA;
            % Calculate xyz of points
            xyz = cloud.XYZ;
            
            % Calculate angle of incidence
            
            internalAngle = zeros(size(xyz,1),1);
            
            for i = 1:length(internalAngle)
                topI = dot(xyz(i,:),SC);
                bottomI = norm(xyz(i,:))*norm(SC);
                
                internalAngle(i) = acos(topI/bottomI);
            end
            
            externalAngle = zeros(size(xyz,1),1);
            for i =1:length(externalAngle)
                topE = dot([SC(1)-xyz(i,1), SC(2)-xyz(i,2), SC(3)-xyz(i,3)],...
                    SC);
                bottomE = norm([SC(1)-xyz(i,1), SC(2)-xyz(i,2), SC(3)-xyz(i,3)])*norm(SC);
                
                externalAngle(i) = acos(topE/bottomE);
            end
            
            AoI = internalAngle + externalAngle;
            
            
            % Calculate variance of alpha and theta
            varTheta = obj.obsStd(2).^2*ones(size(AoI));
            varAlpha = obj.obsStd(3).^2*ones(size(AoI));
            
            varR = ( obj.obsStd(1) .*sec(AoI) ).^2;
            
            varRTA = [varR, varTheta, varAlpha];
            
            varXYZ = zeros(size(varRTA));
            
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
                
                C_x = a* C_L *a';
                
                varXYZ(i,:) = diag(C_x)';
            end
            
        end
        
        function AllClouds = BulkScanPlanes(obj, name)
            % obj.BulkScanPlanes(name) or obj.BulkScanPlanes(planes)
            % this function sequentially scans a set of planes, or reads a
            % file containing the same information. The information of the
            % planes needs to be stored in a 3m x 3 matrix, each column
            % corresponding to x, y and z coordinates of the vertexes, and
            % each set of three rows corresponding to the 3 vertexes of
            % each plane. Only triangular planes are valid, and the normal
            % direction is determined using the cross product of p1->p2 x
            % p1->p3. The points and their variances in both xyz and rho
            % theta alpha are stored in a Cloud object.
            
            AllPlanes = Plane(eye(3));
            if ischar(name)
                var = load(name,'ascii');
            else
                var = name;
            end
            
            if ~isa(var,'Plane')
                
                numPlanes = size(var,1)/3;
                
                for i = 1:numPlanes
                    AllPlanes(i,1) = Plane(var(((3*(i - 1))+1):(3*(i - 1))+3,:));
                end
            else
                AllPlanes = var;
            end
            
            numPlanes = size(AllPlanes,1);
            
            AllClouds = Cloud(zeros(0,3),obj,'angle');
            tcloud = [];
            tvarRTA = [];
            tvarXYZ = [];
            terrRTA = [];
            terrXYZ = [];
            
            for k = 1:numPlanes
                
                
                tempCloud = ScanPlane(obj, AllPlanes(k));
                tcloud = [tcloud; tempCloud.XYZ;];
                tvarRTA = [tvarRTA; tempCloud.varRTA];
                tvarXYZ = [tvarXYZ; tempCloud.varXYZ];
                
                terrRTA = [terrRTA; tempCloud.errRTA];
                terrXYZ = [terrXYZ; tempCloud.errXYZ];

                clear tempCloud
            end
            
            
            
            newTempCloud = Cloud(tcloud,obj,'xyz');
            newTempCloud.varRTA = tvarRTA;
            newTempCloud.varXYZ = tvarXYZ;
            
            newTempCloud.errRTA = terrRTA;
            newTempCloud.errXYZ = terrXYZ;
            AllClouds = AllClouds.combineClouds(newTempCloud);
            
        end
        
        function AllClouds = BulkScanSpheres(obj, name)
            % This function takes in either an array of sphere objects, a
            % matrix containing the xyz of the center and the radius (i.e.
            % a 1x4 vector for each sphere), or a name of a file containing
            % the same information. It returns a cloud object which
            % contains all the scans of each sphere provided. In some
            % cases, spheres will be too far away or too small to have any
            % points fall on their surface, resulting in no data on these
            % points.
            
            AllSpheres = Sphere([0,0,0],0.076);
            
            if ischar(name)
                variable = load(name,'ascii');
            else
                variable = name;
            end
            
            numSpheres = size(variable,1);
            
            if ~isa(variable,'Sphere')
                for i = 1:numSpheres
                    AllSpheres(i,1) = Sphere(variable(i,1:3), variable(i,4));
                end
            else
                AllSpheres = variable;
            end
            
            AllClouds = Cloud(zeros(0,3),obj,'angle');
            tcloud = [];
            tvarRTA = [];
            tvarXYZ = [];
            terrRTA = [];
            terrXYZ = [];
            
            for k = 1:numSpheres
                tempCloud = ScanSphere(obj, AllSpheres(k));
                tcloud = [tcloud; tempCloud.XYZ;];
                tvarRTA = [tvarRTA; tempCloud.varRTA];
                tvarXYZ = [tvarXYZ; tempCloud.varXYZ];
                
                terrRTA = [terrRTA; tempCloud.errRTA];
                terrXYZ = [terrXYZ; tempCloud.errXYZ];
                clear tempCloud
            end
            
            newTempCloud = Cloud(tcloud,obj,'xyz');
            newTempCloud.varRTA = tvarRTA;
            newTempCloud.varXYZ = tvarXYZ;
            
            newTempCloud.errRTA = terrRTA;
            newTempCloud.errXYZ = terrXYZ;
            
            AllClouds = AllClouds.combineClouds(newTempCloud);
        end
        
        function clouds = ScanScene(obj, planesRaw, spheresRaw)
            % This function takes takes in the relevent information for
            % planes and spheres, and returns a cloud which contains the
            % information of all of them. It also combines the clouds such
            % that each point in the cloud has a unique angular location
            % which corresponds to the appropriate, closest surface.
            
            % This mess is to enable many different expressions of planes
            % to be handled. If the variable being passed in is a name of a
            % file, then the file is read, and turned into a matrix. If it
            % is already a matrix, it is then converted into a vector of
            % Planes objects.
            if ischar(planesRaw)
                planesRaw = load(planesRaw,'ascii');
            end
            
            AllPlanes = Plane(eye(3));
            
            if ~isa(planesRaw,'Plane')
                %%%% Delete this later on, for testing purposes only %%%%
                %                 planesRaw = planesRaw/5;
                %%%% This concludes the deleting section %%%%
                
                numPlanes = size(planesRaw,1)/3;
                
                for i = 1:numPlanes
                    AllPlanes(i,1) = Plane(planesRaw(((3*(i - 1))+1):(3*(i - 1))+3,:));
                end
            else
                AllPlanes = planesRaw;
            end
            
            % Same goes for spheres.
            
            if ischar(spheresRaw)
                spheresRaw = load(spheresRaw,'ascii');
            end
            
            AllSpheres = Sphere([0,0,0],0.076);
            
            numSpheres = size(spheresRaw,1);
            
            if ~isa(spheresRaw,'Sphere')
                for i = 1:numSpheres
                    AllSpheres(i,1) = Sphere(spheresRaw(i,1:3), spheresRaw(i,4));
                end
            else
                AllSpheres = spheresRaw;
            end
            
            % This is where we actually do the scanning, having loaded in
            % the coordinates of the planes and spheres
            cloud1 = obj.BulkScanPlanes(AllPlanes);
            cloud2 = obj.BulkScanSpheres(AllSpheres);
            
            clouds = cloud1.combineClouds(cloud2);
            
        end
        
        function [params, C_x] = registerScans(obj, clouds, spheresRaw)
            % [params, C_x] = registerScansMulti(obj, clouds, sphereFileName)
            % This function takes in an array of cloud objects (they do not
            % have to be from the same scan, and in fact, should not be).
            % It also takes in a file name of the spheres, an array of
            % sphere objects, or a matrix of sphere information so that the
            % points which correspond to those spheres can be easily
            % excratced. Should the sphere not be present in the scan, or
            % have an insufficient number of points on its surface, it will
            % not be used in the registration process. This is done
            % automatically.
            
            % clouds is a vector of cloud objects. The scan being used to
            % perform the registration is the home scan.
            
            homeScanIndex = nan;
            for i = 1:length(clouds)
                if(obj.location == clouds(i).scan.location & obj.direction == clouds(i).scan.direction)
                    homeScanIndex = i;
                    break
                end
            end
            
            if isnan(homeScanIndex);
                error('Scan not present in cloud array')
            end
            
            if ischar(spheresRaw)
                spheresRaw = load(spheresRaw,'ascii');
            end
            
            AllSpheres = Sphere([0,0,0],0.076);
            
            numSpheres = size(spheresRaw,1);
            
            if ~isa(spheresRaw,'Sphere')
                for i = 1:numSpheres
                    AllSpheres(i,1) = Sphere(spheresRaw(i,1:3), spheresRaw(i,4));
                end
            else
                AllSpheres = spheresRaw;
            end
            
            
            for i = 1:length(clouds)
                if i < homeScanIndex
                    [movingSC{i}, movingSCvar{i}] = clouds(i).estimateSphereTargets(AllSpheres);
                elseif i == homeScanIndex
                    [HomeSC, HomeSCvar] = clouds(i).estimateSphereTargets(AllSpheres);
                elseif i>homeScanIndex
                    [movingSC{i-1}, movingSCvar{i-1}] = clouds(i).estimateSphereTargets(AllSpheres);
                end
            end
            
            for i = 1:length(movingSC)
                
                a = isnan(movingSC{i});
                b = isnan(HomeSC);
                
                if sum(~b(:,1))<size(HomeSC,1)
                    disp('Warning: Home Scan cannot see all targets')
                end
                
                Key = a|b; Key = ~Key(:,1);
                
                if sum(Key) < 3
                    error('Insufficient Matched Targets between scans. Move scans or add targets')
                end
                
                % Transpose so that points are in columns
                
                HomeSCMulti{i} = deleteRowKey(HomeSC, Key)';
                HomeSCvarMulti{i} = deleteRowKey(HomeSCvar, Key)';
                
                movingSC{i} = deleteRowKey(movingSC{i}, Key)';
                movingSCvar{i} = deleteRowKey(movingSCvar{i}, Key)';
                
            end
            
            params = [];
            for i = 1:length(movingSC)
                [R{i},t{i}]= HornsRegistrationColumns(movingSC{i},HomeSCMulti{i});
                params =[params;[t{i}; atan2(-R{i}(3,2),R{i}(3,3)); asin(R{i}(3,1)); atan2(-R{i}(2,1),R{i}(1,1))]];
            end
            
            initial_params = params;
            max_delta = 1;
            iter = 0;
            
            allDeltas = [];
            
            disp('Beginning Registration')
            while(max_delta > 1e-4 && iter<100)
                A = [];
                
                C_L = [];
                misclosure = [];
                
                for i = 1:length(movingSC)
                    subMatrixA = [];
                    
                    for j = 1:length(movingSC{i})
                        X = movingSC{i}(1,j);
                        Y = movingSC{i}(2,j);
                        Z = movingSC{i}(3,j);
                        
                        w = params(4+(i-1)*6);
                        ph = params(5+(i-1)*6);
                        k = params(6+(i-1)*6);
                        
                        temp = eye(3,6);
                        
                        temp(1,4) = (Z*(cos(w)*sin(k) + cos(k)*sin(w)*sin(ph))) - (Y*(sin(k)*sin(w) - cos(k)*cos(w)*sin(ph)));
                        temp(1,5) = (Y*(cos(k)*cos(ph)*sin(w))) - (Z*(cos(k)*cos(w)*cos(ph)))- (X*(cos(k)*sin(ph)));
                        temp(1,6) = (Y*(cos(k)*cos(w)-sin(k)*sin(w)*sin(ph))) + (Z*(cos(k)*sin(w)+cos(w)*sin(k)*sin(ph))) - (X*(cos(ph)*sin(k)));
                        
                        temp(2,4) = (Z*(cos(k)*cos(w) -sin(k)*sin(w)*sin(ph))) -(Y*(cos(k)*sin(w) + cos(w)*sin(k)*sin(ph)));
                        temp(2,5) = (X*(sin(k)*sin(ph))) + (Z*(cos(w)*cos(ph)*sin(k))) + (Y*(cos(ph)*sin(k)*sin(w)));
                        temp(2,6) = -(Y*(cos(w)*sin(k) + cos(k)*sin(w)*sin(ph))) - (Z*(sin(k)*sin(w)-cos(k)*cos(w)*sin(ph))) - (X*(cos(k)*cos(ph)));
                        
                        temp(3,4) = -(Y*(cos(w)*cos(ph))) - (Z*(cos(ph)*sin(w)));
                        temp(3,5) = (X*cos(ph)) - (Z*(cos(w)*sin(ph))) + (Y*(sin(w)*sin(ph)));
                        
                        subMatrixA = [subMatrixA;temp];
                        
                        clear temp
                    end
                    A = blkdiag(A,subMatrixA);
                end
                
                
                level_constraints = zeros(6*length(movingSC));
                level_con_vector = zeros(6*length(movingSC),1);
                
                for i = 1:length(movingSC)
                    R{i} = (rotz(params(6+((i-1)*6)))*roty(params(5+((i-1)*6)))*rotx(params(4+((i-1)*6))));
                    t{i} = params(1+((i-1)*6):3+((i-1)*6));
                    
                    approxHomeSC = R{i}*movingSC{i} + repmat(t{i},1,size(movingSC{i},2));
                    
                    difference =  approxHomeSC - HomeSCMulti{i};
                    misclosure = [misclosure; difference(:)];
                    
                    C_L = blkdiag(C_L, diag(movingSCvar{i}(:)));
                    
                    level_constraints(4+((i-1)*6),4+((i-1)*6)) = (1.8e-4)^-2;
                    level_constraints(5+((i-1)*6),5+((i-1)*6)) = (1.8e-4)^-2;
                    
                    level_con_vector(4+((i-1)*6)) = params(4+((i-1)*6))/(1.8e-4)^-2;
                    level_con_vector(5+((i-1)*6)) =  params(5+((i-1)*6))/(1.8e-4)^-2;
                    
                end
                
                N = A'*C_L^-1*A + level_constraints;
                u = A'*C_L^-1*misclosure + level_con_vector;
                
                delta_x = -(N\u);
                
                allDeltas = [allDeltas, delta_x];
                
                max_delta = max(abs(delta_x));
                
                params = params + delta_x;
                iter = iter + 1;
                disp(['Iteration Number: ', num2str(iter)])
                if max_delta < 1e-4
                    disp('Registration converged')
                end
            end
            
            C_x = N^-1;
            
            
            
        end
        
        function updatedScans = updateScans(obj, scans, params, C_x, homeScanIndex)
            updatedScans = scans;
            for i = 1:length(scans)
                if i ~= homeScanIndex
                    if i < homeScanIndex
                        updatedScans(i).regParams = params(1+((i-1)*6):6+((i-1)*6));
                        updatedScans(i).regParamVarCovar = C_x(1+((i-1)*6):6+((i-1)*6),1+((i-1)*6):6+((i-1)*6));
                    elseif i > homeScanIndex
                        j = i-1;
                        updatedScans(i).regParams = params(1+((j-1)*6):6+((j-1)*6));                        
                        updatedScans(i).regParamVarCovar = C_x(1+((j-1)*6):6+((j-1)*6),1+((j-1)*6):6+((j-1)*6));
                    end
                    
                end
            end
            
        end
        
    end
    
    
    methods (Static)
        function unwrapped = wrapToPi(ang)
            part = ang + (2*pi * (ang < -pi));
            unwrapped = part - (2*pi * (part>pi) );
        end
        
        
    end
end
