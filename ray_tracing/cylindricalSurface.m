classdef cylindricalSurface < opticalElement_Surface
% cylinder: 
%         R = radius of curvature
%         centAxis = central axis of cylinder (unit vector)
%
%         Notes: The radius should be positive (negative) if the surface is
%         concave (convex) relative to Orientation(:,1).
    properties
        R
        centAxis
    end
    
    methods
        function obj = set.R(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','nonzero'})
            obj.R = value;
        end
        function obj = set.centAxis(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.centAxis = value(:)/norm(value);
        end
        function [d,n] = goToSurf(obj,rayData)
            d = obj.R;
            n = obj.centAxis';
        end
    end
    
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
