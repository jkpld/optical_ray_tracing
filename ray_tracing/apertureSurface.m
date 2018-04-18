classdef apertureSurface < opticalElement_Surface
% apertureSurface  A surface that only allows rays in a center region to
% pass
%
% apertureSurface properties:
%   Shape - Shape of the region of acceptance {circle, rectangle}
%
%   Center - center point of acceptance region
%
%   normVec - normal of plane the aperture is defined on
%
%   Extent - If Shape is circle, then numel(Extent) should = 1 and the
%   value should be the radius where all rays inside the circle may pass
%   through. If Shape is rectangle, then if numel(Extent) = 1, the
%   acceptance region will be a square with edge length Extent(1); if
%   numel(Extent) = 2, the edge length of the acceptance region along
%   Orientation(:,2) will be Extent(1) and the edge length along
%   Orientation(:,3) will be Extent(2).
%
% See also OPTICALELEMENT_SURFACE SPHERICALSURFACE ASPHERICALSURFACE
% PLANESURFACE CYLINDRICALSURFACE CONESURFACE PARABOLICSURFACE

    properties
        Shape = 'circle';
        Center = [0;0;0];
        normVec = [1;0;0];
        Extent = 12.7;
    end

    methods
        function obj = set.Shape(obj, value)
            value = validatestring(value,{'circle','rectangle'});
            obj.Shape = value;
        end
        function obj = set.Center(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Center = value(:);
        end
        function obj = set.normVec(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.normVec = value(:)/norm(value);
        end
        function obj = set.Extent(obj, value)
            validateattributes(value,{'double'},{'finite'})
            numEl = numel(value);
            if numEl>2 || numEl<1
                error('Expected input to be an array with number of elements < 3.');
            end
            obj.Extent = value(:);
        end
        function obj = flipSurf(obj)
        end
        function [d,rayDat,surfNorms, inside] = goToSurf(obj,rayDat,n,~)
            % goToSurf : apertureSurface  propogate rays to an aperture
            % surface and let rays inside the aperture pass through.
            %
            % This function will calculate the distance along a ray that
            % the ray needs to travel to intersect with the plane the
            % aperture is defined on. If the ray is inside of the aperture,
            % it will be allowed to pass; otherwise, it will be blocked (by
            % setting the position of the ray to NaN.
            %
            % surfNorms is not used in this function
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


            d = - (bsxfun(@minus, rayDat(:,1:3),obj.Center')*obj.normVec) ...
                ./ (rayDat(:,4:6)*obj.normVec);

            inside = d > 0;
            d(d<0) = NaN;

            if nargout > 1



                if any(inside)

                    rayDat(inside,1:3) = rayDat(inside,1:3) + bsxfun(@times, rayDat(inside,4:6),d(inside));

                    switch obj.Shape
                        case 'circle'
                            inside = inside & sum(bsxfun(@minus, rayDat(:,1:3),obj.Center').^2,2) < obj.Extent^2;
                        case 'rectangle'
                            if numel(obj.Extent) == 1
                                inside = inside & abs(rayDat(:,2)-obj.Center(2)) < obj.Extent/2 & abs(rayDat(:,3)-obj.Center(3)) < obj.Extent/2;
                            else
                                inside = inside & abs(rayDat(:,2)-obj.Center(2)) < obj.Extent(1)/2 & abs(rayDat(:,3)-obj.Center(3)) < obj.Extent(2)/2;
                            end
                    end
                    if any(inside)
                        rayDat(inside,8) = d(inside).*n(inside);
                    end
                end

                rayDat(~inside,1:3) = NaN;

                if nargout > 2
                    surfNorms = [];
                end
            end

        end

        function createSurfPatch(surfObj, opticalEl, surfNum, ~)
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

            % Get cross section boundry points
            [~, ~, Yb, Zb, ~, ~] = opticalElement_Surface.generateSurfMesh(opticalEl);

            uv = cross([1;0;0],surfObj.normVec);

            ang = asin(norm(uv));
            uv = uv/norm(uv);

            rot = eye(4);
            if ang~=0
                rot = makehgtform('axisrotate',uv,ang);
            end

            cent = surfObj.Center;

            cent = rot*[cent;1];

            Yb = Yb - cent(2);
            Zb = Zb - cent(3);

            N = surfObj.numSurfPoints;
            switch surfObj.Shape
                case 'circle'
                    r = surfObj.Extent;
                    th = linspace(0,2*pi,N);

                    Yi = r*cos(th);
                    Zi = r*sin(th);

                case 'rectangle'
                    l = surfObj.Extent;
                    if numel(l) == 1
                        Yi = [1;-1;-1;1;1]*l;
                        Zi = [1;1;-1;-1;1]*l;
                    else
                        Yi = [1;-1;-1;1;1]*l(1);
                        Zi = [1;1;-1;-1;1]*l(2);
                    end


            end

            [Yb,Zb] = poly2cw(Yb,Zb);
            Y = [Yb(:);NaN;Yi(:)];
            Z = [Zb(:);NaN;Zi(:)];

            [f,v] = poly2fv(Y,Z);
            % Create 4x4 rotation/translation matrix to position the
            % element to have normal vector normVec and aperature center
            % Center.


            rot(1:3,4) = surfObj.Center;

%             X = [zeros(length(Yb),1);NaN;zeros(length(Yi),1)];
%             Y = [Yb(:);NaN;Yi(:)];
%             Z = [Zb(:);NaN;Zi(:)];

            XYZ = rot*[zeros(1,size(v,1));v';ones(1,size(v,1))];

            v = XYZ(1:3,:)';


            patch('Vertices', v, ...
                  'Faces', f, ...
                  'FaceColor', [0.2 0.2 0.2], ...
                  'Parent', group, ...
                  'FaceAlpha', 0.9, ...
                  'FaceLighting','none',...
                  'EdgeColor', 'none', ...
                  'LegendDisplay', 'off', ...
                  'HitTest', 'off',...
                  'Tag','apertureSurface');

            opticalEl.Surfaces(surfNum).boundingEdge = XYZ(1:3,:)'; % Here we cannot just save it to surfObj, because surfObj is not a handle and we are not returning surfObj. So we can just change the actual surfObj inside the handle class opticalEl.
        end
    end

end

%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
