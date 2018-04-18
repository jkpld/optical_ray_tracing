classdef opticalElement_Surface < matlab.mixin.Heterogeneous
    % Parent class for the different surface classes. The surface classes are
    %   sphericalSurface
    %   asphereicalSurface
    %   cylindricalSurface
    %   planeSurface
    %   conicalSurface
    %   aperatureSurface
    %   nullSurface

    properties(Hidden = true)
        clearAperture    = []; % Clear aperature of the surface in mm

        % opticalCoating % This is not implemented, but if it is,
        % opticalCoating should probably be a class object
    end

    properties(Constant, Hidden = true)
        numSurfPoints = 51;
        baseFaceColor = [.816, .902, .996]*1.3;
%         plotOptions = struct('FaceAlpha', 0.4, ...
%                          'FaceLighting','phong',...
%                          'BackFaceLighting', 'lit', ...
%                          'EdgeColor', 'none', ...
%                          'LegendDisplay', 'off', ...
%                          'HitTest', 'off')
    end
    
    properties(Hidden = true)
        boundingEdge = []; % n x 3 array giving the bounding edge [x, y, z]
    end
    
    methods (Abstract)
        goToSurf(obj,rayData,n,opticalEl)
        createSurfPatch(obj,opticalEl,surfNum)
        flipSurf(obj)
    end
    
    methods (Static, Sealed, Access = protected)
        function default_object = getDefaultScalarElement
            default_object = nullSurface;
        end
        
        function [Yp, Zp, Yb, Zb, Yc, Zc] = generateSurfMesh(opticalEl)
            N = opticalEl.Surfaces(1).numSurfPoints;
            switch opticalEl.crossSection
                case 'circle'
                    r = opticalEl.crossSectionSize(1);
                    th = linspace(0,pi/2,ceil(N/4));
                    th = [th(1:end-1), th(1:end-1) + pi/2, th(1:end-1) + pi, th + 3*pi/2];
                    [Rp,Thp] = meshgrid(linspace(0,r,ceil(N/2)),th);
                    Yp = Rp.*cos(Thp);
                    Zp = Rp.*sin(Thp);
                    
                    Yb = r*cos(th'); % y points for boarder edge
                    Zb = r*sin(th'); % z points for boarder edge
                    
                    Yc = [linspace(-r,r,N)'; NaN; zeros(N,1)];
                    Zc = flip(Yc,1);
                    
                case 'rectangle'
                    if numel(opticalEl.crossSectionSize) == 1
                        l = opticalEl.crossSectionSize(1);
                        l = linspace(-l/2,l/2,N);
                        [Yp,Zp] = meshgrid(l,l);
                    else % numel = 2
                        l = opticalEl.crossSectionSize;
                        [Yp,Zp] = meshgrid(linspace(-l(1)/2,l(1)/2,N),linspace(-l(2)/2,l(2)/2,N));
                    end

                    Yb = [Yp(1,:)';Yp(2:end,end);Yp(end,(end-1):-1:1)';Yp((end-1):-1:1,1)];
                    Zb = [Zp(1,:)';Zp(2:end,end);Zp(end,(end-1):-1:1)';Zp((end-1):-1:1,1)];
                    
                    Yc = [linspace(-l(1)/2,l(1)/2,N)'; NaN; zeros(N,1)];
                    Zc = [zeros(N,1); NaN; linspace(-l(2)/2,l(2)/2,N)'];
            end
        end
    end
    
    methods (Sealed)
        function out = isa(obj,classType)
            out = false(size(obj));
            for i = 1:numel(obj)
                out(i) = builtin('isa',obj(i),classType);
            end
        end
        
        function flipped = flip(obj)
            % flip  Flip the surface array and change the convexivity

            flipped = builtin('flip',obj);

            for i = 1:numel(flipped)
                flipped(i) = flipSurf(flipped(i));
            end
        end
        
    end
    
    methods (Static, Sealed) % To use a static method you must call opticalElement_Surface.(methodName)
        
        function options = displayOptions(obj, material)

            options = struct('FaceAlpha', 0.4, ...
                         'FaceLighting','gouraud',...
                         'BackFaceLighting', 'reverselit', ...
                         'EdgeColor', 'none', ...
                         'LegendDisplay', 'off', ...
                         'HitTest', 'off');
            props = [0.3 0.3 1 25 0.5];
            switch material
                case {'Aluminum','Al'}
                    options.FaceColor = [173 178 189]/255;
                    options.FaceAlpha = 0.9;
                    options.AmbientStrength = props(1);
                    options.DiffuseStrength = props(2);
                    options.SpecularStrength = props(3);
                    options.SpecularExponent = props(4);
                    options.SpecularColorReflectance = props(5);
                case {'Silver','Ag'}
                    options.FaceColor = [192 192 195]/255;
                    options.FaceAlpha = 0.9;
                    options.AmbientStrength = props(1);
                    options.DiffuseStrength = props(2);
                    options.SpecularStrength = props(3);
                    options.SpecularExponent = props(4);
                    options.SpecularColorReflectance = props(5);
                case {'Gold','Au'}
                    options.FaceColor = [212 175 55]/255;
                    options.FaceAlpha = 0.9;
                    options.AmbientStrength = props(1);
                    options.DiffuseStrength = props(2);
                    options.SpecularStrength = props(3);
                    options.SpecularExponent = props(4);
                    options.SpecularColorReflectance = props(5);
                otherwise
                    options.FaceColor = obj.baseFaceColor/refractiveIndex(material,532);
            end

        end
        
        function [vnew,fnew] = removeDubs(V,F)
            % From comment here
            % http://www.mathworks.com/matlabcentral/fileexchange/29986-patch-slim--patchslim-m

            [vnew, ~, indexn] = unique(V, 'rows'); 
            fnew = indexn(F);

            %remove nonsens faces 
            numfaces = (1:size(fnew,1))';

            e1=fnew(:,1)-fnew(:,2); 
            e2=fnew(:,1)-fnew(:,3); 
            e3=fnew(:,2)-fnew(:,3);

            e1=[e1 numfaces]; 
            e2=[e2 numfaces]; 
            e3=[e3 numfaces];

            e1=e1(e1(:,1)==0,2); 
            e2=e2(e2(:,1)==0,2); 
            e3=e3(e3(:,1)==0,2);

            fnew(vertcat(e1,e2,e3),:)=[];
        end
        
        
    end
end
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
