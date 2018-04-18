classdef sphericalSurface < opticalElement_Surface
    % sphere: parameters, size(parameters) = [4, 1]
    %         parameters(1) = R, R~=0 && abs(R)~=inf
    %         parameters(2:4) = center of sphere
    %
    %         Notes: (1) The radius should be positive
    %         (negative) if the surface is convex (concave)
    %         relative to Orientation(:,1). (2) The center of
    %         the sphere is relative to the elements Center.
    properties
        R
        Center
    end
    
    methods
        function obj = set.R(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','nonzero'})
            obj.R = value;
        end
        function obj = set.Center(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Center = value(:);
        end
        
        function [d,rayDat,surfNorms,inside] = goToSurf(obj, rayDat, n, el)
            % goToSurf : sphericalSurface  propogate rays to sphere and
            % calculate the sphere normals at the ray intersection points
            %
            % This function will calculate the distance along a ray that
            % the ray needs to travel to intersect with a spherical
            % surface. The rays are then propogated to the intersection
            % points, and the surface normals of the spherical surface at
            % the intersection points calculated. The output will be the
            % distances, the new array of rayDat, and the spherical surface
            % normals.
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
            % surfNorms - surface normals as the points of intersection
            %
            % inside - logical array of size [size(rayDat,1),1]. Elements
            % are true if the ray intersects the surface and the point of
            % intersection is inside of the elements crossSection.
            
            % Written by James Kapaldo, 2016, 01, 22
            % 
            % Modification History:
            % 2016, 01, 23 : added distance to the output and some nargout
            % logic
            % 2016, 01, 24 : added the optical element as an input, added
            % 'inside' as an output, and now calculate if a ray that
            % inersects the surface also is inside of the elements
            % crossSection.

            radius = obj.R;
            sphereCent = obj.Center;
            y = bsxfun(@minus, rayDat(:,1:3),sphereCent');

            y_dot_uv = sum(y.*rayDat(:,4:6),2);
            
            underSqrt = y_dot_uv.^2 - (sum(y.^2,2) - radius^2);
            
            good = underSqrt>=0;

            d = NaN(size(rayDat,1),1);
            
            d(good) = -y_dot_uv(good) - sign(radius)*sqrt(underSqrt(good));
            d(d<0) = NaN;
            
            if nargout > 1
                rayDat(:,8) = d.*n;
                rayDat(:,1:3) = rayDat(:,1:3) + bsxfun(@times, rayDat(:,4:6), d);
                
                [rayDat, inside] = removeRaysOutsideCrossSection(el,rayDat);
                
                if nargout > 2
                    surfNorms(length(inside),3) = 0;
                    
                    surfNormsIn = bsxfun(@minus, rayDat(inside,1:3), sphereCent');
                    surfNormsIn = bsxfun(@rdivide, surfNormsIn, sqrt(sum(surfNormsIn.*surfNormsIn,2)));
                    
                    surfNorms(inside,:) = surfNormsIn;
                end
            end

        end
        
        function obj = flipSurf(obj)
            obj.R = -obj.R;
%             obj.Center = -obj.Center;
        end
        
        function createSurfPatch(surfObj, opticalEl, surfNum, material)
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
            
            Xp = sign(surfObj.R)*sqrt(surfObj.R^2 - ((Yp-surfObj.Center(2)).^2 + (Zp-surfObj.Center(3)).^2)) + surfObj.Center(1);
            Xb = sign(surfObj.R)*sqrt(surfObj.R^2 - ((Yb-surfObj.Center(2)).^2 + (Zb-surfObj.Center(3)).^2)) + surfObj.Center(1);
            Xc = sign(surfObj.R)*sqrt(surfObj.R^2 - ((Yc-surfObj.Center(2)).^2 + (Zc-surfObj.Center(3)).^2)) + surfObj.Center(1);
            
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
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
