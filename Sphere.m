classdef Sphere
    properties
        center
        radius
    end
    
    methods
        function obj = Sphere(x,r)
            obj.center = x;
            obj.radius = r;
        end
    end
end