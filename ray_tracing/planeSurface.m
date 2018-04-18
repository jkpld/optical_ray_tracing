classdef planeSurface < opticalElement_Surface
% plane:
%       Point = a point on the plane
%       normVec = unit normal vector of the plane
    properties
        Point
        normVec
    end
    
    methods
        function obj = set.Point(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Point = value(:);
        end
        function obj = set.normVec(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            if value(1) == 0
                error('planeSurface:badNormVec','The normal vector may not be perpendicular to the x-axis.')
            end
            obj.normVec = value(:)/norm(value);
        end
        function obj = flipSurf(obj)
        end
        function [d,rayDat,surfNorms, inside] = goToSurf(obj,rayDat,n,el)
            % goToSurf : planeSurface  Propogate rays to a plane surface
            % and get the surface normals
            %
            % This function will calculate the distance along a ray that
            % the ray needs to travel to intersect with a plane
            % surface. The rays are then propogated to the intersection
            % points, and the surface normals of the spherical surface at
            % the intersection points calculated. The output will be the
            % distances, the new array of rayDat, and the spherical surface
            % normals, and an array inside that is true for all rays inside
            % the element's crossSection
            %
            %
            % Inputs:
            % obj - A surface object
            % rayDat - n x 8 matrix of ray data, where n is the number
            %          of rays. The 10 columns are
            %          x,y,z,uv_x,uv_y,uv_z,w,l
            %
            %          x,y,z -- starting position of ray
            %          uv_x,uv_y,uv_z -- unit vector giving the 
            %                            direction of the ray
            %          w -- wavelength of ray
            %          l -- total optical path length the ray has 
            %               traveled
            %
            % n - refractive index of the medium the rays are traveling in
            %
            % el - The opticalElement3D that the surface belongs to
            %
            % Outputs:
            % d - The distance needed to travel along the rays.
            % rayDat - The ray data at the points of intersection
            %
            % surfNorms - not used. Empty will be returned.
            %
            % inside - logical array of size [size(rayDat,1),1]. Elements
            % are true if the ray passes through the aperture and false
            % otherwise.
            
            % Written by James Kapaldo, 2016, 01, 24
            % 
            % Modification History:
            % 
            
            % n.(r-r0) = 0  , n = normVec, r0 = Center, r = equation of
            % line
            % normVec.(a+d*uv - center) = 0
            % d = normVec(center - a)/normVec.uv
            
            
            d = - (bsxfun(@minus, rayDat(:,1:3),obj.Point')*obj.normVec) ...
                ./ (rayDat(:,4:6)*obj.normVec);
            
            d(d<0) = NaN;
            
            if nargout > 1
                rayDat(:,8) = d.*n;
                rayDat(:,1:3) = rayDat(:,1:3) + bsxfun(@times, rayDat(:,4:6), d);
                
                [rayDat, inside] = removeRaysOutsideCrossSection(el,rayDat);
                
                if nargout > 2
                    if any(inside)
                        surfNorms(length(inside),3) = 0;                    
                        surfNorms(inside,:) = repmat(obj.normVec',sum(inside),1);
                    else
                        surfNorms = [];
                    end
                end
            end
        end
        
        function createSurfPatch(surfObj, opticalEl, surfNum,material)
            % Assume opticalEl is a valid element (no messed up surfaces)
            % Assume opticalEl.Center = 0 and opticalEl.Orientation =
            % eye(3). Then after creating the entire optical element, we
            % will just transform the group
            
            props = fieldnames(surfObj)';
            for i = 1:numel(props)
                if isempty(surfObj.(props{i}))
                    error('opticalElement_Surface:undefinedParameter','A parameter of %s is undefined and is, therefore, not valid.',class(surfObj))
                end
            end
            
            % Create a group to hold the different surface plot objects --
            % the primary patach, the 3 lines to help visualization, ...
            group = matlab.graphics.primitive.Group;
            group.Parent = opticalEl.ShapeGroupPrimary;
            
            % create grid of cross section points to evaluate the surface
            % on
            [Yp, Zp, Yb, Zb, Yc, Zc] = opticalElement_Surface.generateSurfMesh(opticalEl);

            % Now project the surface points, boundary points, and cross
            % points to the plane.
            % Equation of plane : n.(r-r0)=0
            n = surfObj.normVec;
            r0 = surfObj.Point;
            
            Xp = -(n(2)*(Yp-r0(2))+n(3)*(Zp-r0(3)))/n(1) + r0(1);
            Xb = -(n(2)*(Yb-r0(2))+n(3)*(Zb-r0(3)))/n(1) + r0(1);
            Xc = -(n(2)*(Yc-r0(2))+n(3)*(Zc-r0(3)))/n(1) + r0(1);
            
            [f,v] = surf2patch(Xp,Yp,Zp);
            [v,f] = opticalElement_Surface.removeDubs(v,f);
            
            options = surfObj.displayOptions(surfObj, material);
            
            patch('Vertices', v, ...
                  'Faces', f, ...
                  'Parent', group, ...
                  options);

            line(Xb,Yb,Zb,'Color',[0.2 0.2 0.2],'Parent',group,'Tag','boundryEdge')
            line(Xc,Yc,Zc,'Color',[0.2 0.2 0.2],'Parent',group,'Tag','cross')
            
            opticalEl.Surfaces(surfNum).boundingEdge = [Xb(:),Yb(:),Zb(:)]; % Here we cannot just save it to surfObj, because surfObj is not a handle and we are not returning surfObj. So we can just change the actual surfObj inside the handle class opticalEl.
        end

    end
    
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
