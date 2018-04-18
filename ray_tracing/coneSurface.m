classdef coneSurface < opticalElement_Surface
% conical surface: (think axicon)
%       Angle = angle of the cone
%       normVec = unit normal vector of the plane (the vector the cone is circularly symetric about
%       Vertex = the cone vertex location
    properties
        Angle
        Vertex
        normVec
    end
    
    methods
        function obj = set.Angle(obj, value)
            validateattributes(value,{'double'},{'finite','numel',1,'<',pi,'>',-pi})
            obj.Angle = value(:);
        end
        function obj = set.Vertex(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Center = value(:);
        end
        function obj = set.normVec(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.normVec = value(:)/norm(value);
        end
        

        function [d,n] = goToSurf(obj,rayData)
            d = obj.Vertex';
            n = obj.Angle;
        end
    end
    
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
