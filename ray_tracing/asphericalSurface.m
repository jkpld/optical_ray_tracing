classdef asphericalSurface < opticalElement_Surface
% asphericalSurface  An aspherical surface
%
% asphericalSurface properties:
%   Center - center of the asphere (the apex of the asphere)
%
%   R - Base radius of the asphere. If R > 0, it is concave; if R < 0, it
%   is convex. (finite, nonzero, scalar)
%
%   k - Conic section (finite, scalar)
%
%   A - Deviation coefficients
%
% See also OPTICALELEMENT_SURFACE SPHERICALSURFACE PARABOLICSURFACE
% PLANESURFACE CYLINDRICALSURFACE CONESURFACE APERTURESURFACE

    properties
        Center
        R
        k
        A
    end
    
    methods
        function obj = set.Center(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Center = value(:);
        end
        function obj = set.R(obj, value)
            validateattributes(value,{'double'},{'finite','scalar','nonzero'})
            obj.R = value;
        end
        function obj = set.k(obj, value)
            validateattributes(value,{'double'},{'finite','scalar'})
            obj.k = value;
        end
        function obj = set.A(obj, value)
            validateattributes(value,{'double'},{'finite','vector'})
            idx = find(value,1,'last');
            if idx~=length(value)
                value(idx+1:end) = [];
            end
            obj.A = value(:);
        end
        function obj = flipSurf(obj)
%             obj.R = -obj.R;
        end
        function [d, rayDat, surfNorms, inside] = goToSurf(obj,rayDat,n,el)
            % goToSurf : asphericSurface  propogate rays to an aspheric
            % surface and get the surface normals at the intersection
            % points.
            %
            % This function will calculate the distance along a ray that
            % the ray needs to travel to intersect with an aspheric surface
            % (d). The rays are then propogated to the intersection points,
            % and the surface normals of the aspheric surface at the
            % intersection points calculated. The output will be the
            % distances, the new array of rayData, the aspherical
            % surface normals, and an array inside that is true for the
            % rays that are inside of the elements cross section.
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
            
            % Written by James Kapaldo, 2016, 01, 25
            % 
            % Modification History:
            
            
            function [fun, dfun] = F(d)
                oneToA = 1:length(obj.A);
                signR = sign(obj.R);
                absR = abs(obj.R);
                r = sqrt((rayDat(:,2)+d.*rayDat(:,5)).^2 + (rayDat(:,3)+d.*rayDat(:,6)).^2);
                
                gama = sqrt(1-(1+obj.k)*r.^2/absR^2);

                f = -signR * ( r.^2./(absR .* (1+gama)) + bsxfun(@power,r,2*oneToA)*obj.A ) + obj.Center(1) - rayDat(:,1) - d.*rayDat(:,4);
                
                if nargout > 1
                    dfdr = -signR * ( 2./(absR.*(1+gama)) + r.^2.*(1+obj.k)./(absR^3 .* gama .* (1+gama).^2) + bsxfun(@power,r,2*oneToA - 2)*(2*obj.A.*oneToA') );
                    drdd = (d.*(rayDat(:,5).^2 + rayDat(:,6).^2) + rayDat(:,2).*rayDat(:,5) + rayDat(:,3).*rayDat(:,6));
%                     [r,f,dfdr,drdd]
                    fun = sum(f.^2);
                    dfun = 2*f.*(dfdr.*drdd-rayDat(:,4));
                else
                    fun = f;
                end
            end
            
            % Create an intial guess for the distance needed to travel.
            numRys = size(rayDat,1);
            
            if ~isempty(obj.boundingEdge)
                x_trial = mean([mean(obj.boundingEdge(:,1)), obj.Center(1)]);
            else
                x_trial = obj.Center(1);
            end
            
            d_trial = (x_trial - rayDat(:,1)) ./ (rayDat(:,4));

            % Two possible methods of solving:
            
            % method 1: create the objective function for a single ray and
            % then inside a forloop solve for the distance to the
            % aspherical surface for each ray
            
%             d(numel(numRys),1) = 0;
            inside = true(numRys,1);
%             oR = obj.R;
%             ok = obj.k;
%             oA = obj.A;
%             oC1 = obj.Center(1);
            
%             for i = 1:numel(d_trial)
% %                 F = @(d) asphere(rayDat(i,2)+d*rayDat(i,5), rayDat(i,3)+d*rayDat(i,6), obj.R, obj.k, obj.A) + obj.Center(1) - rayDat(i,1) - d*rayDat(i,4);
%                 try
%                     d(i) = fzero(@(d) F(d, rayDat(i,1:6),oR,ok,oA,oC1), d_trial(i));
%                 catch
%                     d(i) = NaN;
%                     inside(i) = false;
%                 end
%             end
            
            
            
            % method 2: create an objective function to minimize all rays
            % at the same time using fminunc. Use fminunc because I will
            % not use any contraints, and I know what the gradient is.
            % Also, put in that the hessian is diagonal since each ray is
            % independent.

            options = optimoptions(@fminunc,'Algorithm','trust-region','GradObj','on','HessPattern',speye(numRys),'TolX',1e-14,'TolFun',1e-14,'Display','off');
            d = fminunc(@F, d_trial,options);

            inside(d<0) = false;
            d(~inside) = NaN;
            
            d = d(:);
            
            if nargout > 1
                rayDat(inside,8) = d(inside).*n(inside);
                rayDat(inside,1:3) = rayDat(inside,1:3) + bsxfun(@times, rayDat(inside,4:6), d(inside));
                
                [rayDat(inside,:,:), in] = removeRaysOutsideCrossSection(el,rayDat(inside,:,:));
                inside(inside) = in;
                
                if nargout > 2
                    
                    oneToA = 1:length(obj.A);

                    r = sqrt(rayDat(inside,2).^2 + rayDat(inside,3).^2);

                    gama = sqrt(1-(1+obj.k)*r.^2/obj.R^2);

                    dfdr = -sign(obj.R) * ( 2./(abs(obj.R).*(1+gama)) + r.^2.*(1+obj.k)./(abs(obj.R)^3 .* gama .* (1+gama).^2) + bsxfun(@power,r,2*oneToA - 2)*(2*obj.A.*oneToA') );
                    drdy = rayDat(inside,2);
                    drdz = rayDat(inside,3);

                    % Get surface normals
                    surfNorms(numel(inside),3) = 0;

                    surfNormsIn = [-ones(numel(r),1), dfdr.*drdy, dfdr.*drdz];
                    surfNormsIn = bsxfun(@rdivide, surfNormsIn, sqrt(sum(surfNormsIn.*surfNormsIn,2)));
                    
                    surfNorms(inside,:) = surfNormsIn;
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

            Xp = asphere(Yp-surfObj.Center(2),Zp-surfObj.Center(3),surfObj.R,surfObj.k,surfObj.A) + surfObj.Center(1);
            Xb = asphere(Yb-surfObj.Center(2),Zb-surfObj.Center(3),surfObj.R,surfObj.k,surfObj.A) + surfObj.Center(1);
            Xc = asphere(Yc-surfObj.Center(2),Zc-surfObj.Center(3),surfObj.R,surfObj.k,surfObj.A) + surfObj.Center(1);
            
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

function z = asphere(x,y,R,k,A)
% aphsere  Create a 2d aspheric surface giving x and y arrays and the
% aspheric surface parameters: R, the radius; k, the conic section; and A,
% the deviation coefficients.

z = -sign(R)*reshape((x(:).^2 + y(:).^2)./(R*(1+sqrt(1-(1+k)*(x(:).^2 + y(:).^2)/R^2))) + bsxfun(@power,sqrt(x(:).^2 + y(:).^2),2*(1:length(A)))*A(:),size(x));

end


% Define function for intersecting the asphere
% function [fun, dfun] = F(d,rayDat,R,k,A,c1)
% oneToA = 1:length(A);
% signR = sign(R);
% R = abs(R);
% r = sqrt((rayDat(:,2)+d*rayDat(:,5)).^2 + (rayDat(:,3)+d*rayDat(:,6)).^2);
% 
% gama = sqrt(1-(1+k)*r.^2/R^2);
% 
% f = -signR * ( r.^2/(R .* (1+gama)) + bsxfun(@power,r,2*oneToA)*A ) + c1 - rayDat(:,1) - d*rayDat(:,4);
% 
% if nargout > 1
%     dfdr = -signR * ( 2*r./(R.*(1+gama)) + r.^3.*(1+k)./(R^3 .* gama .* (1+gama).^2) + bsxfun(@power,r,2*oneToA - 1)*(2*A.*oneToA') );
%     drdd = (d*(rayDat(:,5).^2 + rayDat(:,6).^2) + rayDat(:,2).*rayDat(:,5) + rayDat(:,3).*rayDat(:,6))./r;
%     
%     fun = sum(f.^2);
%     dfun = 2*f.*(dfdr.*drdd-rayDat(:,4));
% else
%     fun = f;
% end
% end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
