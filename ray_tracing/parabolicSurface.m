classdef parabolicSurface < opticalElement_Surface
% parabolicSurface  A converging reflective parabolic surface 
%
% The collimated direction will be along the optical axis Orientation(:,1)
% and the reflected direction will be in the plane made by Orientation(:,1)
% and Orientation(:,2).
%
% parabolicSurface properties:
%   PFL - Focal length of parent parabola (parent focal length)
%
%   EFL - Effective focal length, or reflection focal length
%
%   Angle - The reflection angle of the parabola (in degrees)
%
%   baseDist - The distance along the optical axis between the focal point
%   and the Center.
%
% See also OPTICALELEMENT_SURFACE SPHERICALSURFACE ASPHERICALSURFACE
% PLANESURFACE CYLINDRICALSURFACE CONESURFACE APERTURESURFACE

    properties
        PFL
        EFL
        Angle
        baseDist
    end
    
    methods
        function obj = set.PFL(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','nonzero'})
            obj.PFL = value(:);
        end
        function obj = set.EFL(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','nonzero'})
            obj.EFL = value(:);
        end
        function obj = set.Angle(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','>=',0,'<',180})
            obj.Angle = value;
        end
        function obj = set.baseDist(obj, value)
            validateattributes(value,{'double'},{'finite','scalar'})
            obj.baseDist = value(:);
        end
        function obj = flipSurf(obj)
        end
        function [d,rayDat, surfNorms, inside] = goToSurf(obj,rayDat,n,el)
            % goToSurf : parabolicSurface  propogate rays a parabolic
            % surface and get the surface normals at the intersection
            % points.
            %
            % This function will calculate the distance along a ray that
            % the ray needs to travel to intersect with a parabolic surface
            % (d). The rays are then propogated to the intersection points,
            % and the surface normals of the parabolic surface at the
            % intersection points calculated. The output will be the
            % distances, the new array of rayData, and the spherical
            % surface normals.
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
            
            % Written by James Kapaldo, 2016, 01, 24
            % 
            % Modification History:
            
            r0 = [obj.EFL*sind(obj.Angle); 0];

            lx = rayDat(:,1)+obj.PFL;
            lr = rayDat(:,2:3);
            
            ux = rayDat(:,4);
            ur = rayDat(:,5:6);
            
            a = sum(ur.*ur,2);
            b = sum(lr.*ur,2) - ur*r0 - 2*obj.PFL*ux;
            c = sum(lr.*lr,2) + r0'*r0 - 4*obj.PFL*lx - 2*lr*r0;
            
            aeq0 = a == 0;
            
            d(numel(a),1) = 0;
            inside = false(numel(a),1);
            
            if any(aeq0)
                d(aeq0) = -c(aeq0)./(2*b(aeq0));

                int1 = rayDat(aeq0,1:3) + bsxfun(@times, rayDat(aeq0,4:6), d(aeq0));
                [~,inaeq0] = removeRaysOutsideCrossSection(el,int1);

                inside(aeq0) = inaeq0;
            end
            
            if any(~aeq0)
                d_p = (-b(~aeq0) + sqrt(b(~aeq0).^2-a(~aeq0).*c(~aeq0)))./a(~aeq0);
                d_m = (-b(~aeq0) - sqrt(b(~aeq0).^2-a(~aeq0).*c(~aeq0)))./a(~aeq0);

                int1 = rayDat(~aeq0,1:3) + bsxfun(@times, rayDat(~aeq0,4:6), d_p);
                int2 = rayDat(~aeq0,1:3) + bsxfun(@times, rayDat(~aeq0,4:6), d_m);

                [~,in1] = removeRaysOutsideCrossSection(el,int1);
                [~,in2] = removeRaysOutsideCrossSection(el,int2);

                inaneq0 = xor(in1,in2);  % If there is two intersections inside the crossSection, then the ray must be traveling through the back of the mirror

                dq(numel(inaneq0),1) = 0;
                dq(inaneq0 & in1) = d_p(inaneq0 & in1);
                dq(inaneq0 & in2) = d_m(inaneq0 & in2);

                d(~aeq0) = dq;
                inside(~aeq0) = inaneq0;
            end
            
            inside = inside & d >= 0;
            
            d(d<0) = NaN;
                       
            % Get intersection points
            rayDat(inside,1:3) = rayDat(inside,1:3) + bsxfun(@times, rayDat(inside,4:6), d(inside));

            % Get surface normals
            surfNorms(numel(inside),3) = 0;

            surfNormsIn = [ones(sum(inside),1), ...
                -1/(2*obj.PFL) * (rayDat(inside,2)-r0(1)), ...
                -1/(2*obj.PFL) * (rayDat(inside,3)-r0(2))];
            surfNormsIn = bsxfun(@rdivide, surfNormsIn, sqrt(sum(surfNormsIn.*surfNormsIn,2)));

            surfNorms(inside,:) = surfNormsIn;
            
            % Lastly make sure the rays are pointing towards the mirror
            % (and not from the back).
            inside(inside) = inside(inside) & sum(rayDat(inside,4:6).*surfNormsIn,2) < 1;
            
            surfNorms(~inside,:) = 0;
            d(~inside) = NaN;
            
            rayDat(inside,8) = d(inside).*n(inside);
            rayDat(~inside,[1:3,8]) = NaN;

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

            offset = surfObj.EFL*sind(surfObj.Angle);
            
            % x = r^2/(4*p)
            Xp = ((Yp-offset).^2 + Zp.^2 - offset^2)/(4*surfObj.PFL) + surfObj.baseDist;
            Xb = ((Yb-offset).^2 + Zb.^2 - offset^2)/(4*surfObj.PFL) + surfObj.baseDist;
            Xc = ((Yc-offset).^2 + Zc.^2 - offset^2)/(4*surfObj.PFL) + surfObj.baseDist;
            
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
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
