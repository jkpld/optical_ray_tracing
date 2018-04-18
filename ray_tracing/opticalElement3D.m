classdef opticalElement3D < matlab.mixin.Copyable
% opticalElement3D   A 3D optical element for ray tracing 
% An object from this class will store the necessary information for
% creating an optical element and ray tracing that element.
%
% opticalElement3D Properties:
%   Primary properties, 
%       crossSection - The cross sectional shape of the element
%       {'circle','rectangle'}
%
%       crossSectionSize - Cross sectional size of the element. If a
%       'circle', than crossSectionSize(1) gives the radius. If a
%       'rectangle', then, if numel(crossSecionSize) = 1,
%       crossSectionSize(1) gives the edge length for the square; if
%       numel(crossSectionSize) = 2, than crossSectionSize(1) gives the
%       edge length along Orientation(:,2), and crossSectionSize(2) gives
%       the edge length along Orientation(:,3).
%
%       Orientation - The local coordinate system of the element. Must be a
%       3x3 orthogonal matrix. The first column gives the the optical axis
%       of the element. The second and third columns give two
%       perpendicular directions in the crossSection plane. When setting
%       this property, you may give a 3x2 array where the first column is a
%       unit vector giving the optical axis and the second column is a
%       perpendicular unit vector. If it is assigned in this way, then the
%       third direction will be automatically determined from
%       cross(colomn1,column2).
%
%       Origin - A vector giving the location of the origin of the element.
%
%       Surfaces - An array of <a href="matlab: help
%       opticalElement_Surface">opticalElement_Surface</a> describing both
%       the surfaces that make up the optical element and any additional
%       aperatures (pupil, iris). The order in which the surfaces are
%       traced is determined by the sequentialRayTrace property.
%
%       sequentialRayTrace - (logical) If true, then the surfaces of the
%       element will be traced sequentially, either forwards (Surfaces(1)
%       -> Surfaces(2) -> ... -> Surfaces(end)) or backwards (Surfaces(end)
%       -> Surfaces(end-1) -> ... -> Surfaces(1)) depending on if a ray
%       hits Surface(1) or Surface(end) first. If sequentialRayTrace is
%       false, then the element will be traced non-sequentially by
%       determining which surface in the array Surfaces a ray hits first,
%       and then following the ray until it would exit the element's
%       Extent.
%
%       Material - A cell array of strings giving the materials in the
%       element. The number of materials must be equal to the number of
%       non-aperature surfaces in Surfaces minus 1.
%
%       Extent - This is a derived array (3x8) of the 8 vertices giving the
%       bounding volume of the element.
%
%       ExtentAA - This is a derived array (2x3) giving the minimum vertex
%       (first row) and maximum vertex (second row) of the element's axis
%       aligned bounding box.
%
%       Matrix - This is a derived 4x4 matrix that stores the coordinate
%       transformation of the element.
%
%       ShapeGroup - This is a derived Transform group that contains the
%       handles to all of the patches and lines that make up the element.
%       By default the ShapeGroup has no Parent. To display the element, se
%       the ShapeGroup.Parent to an axis handle.
%       
%   Additional meta-data properties,
%       Type - An identification type for the optical element
%
%       SubType - An identification sub-type for the optical element
%
%       Focal - Focal length of the element. This will normally be a 2
%       element vector such that the first element is the focal length and
%       the second is the back focal length.
%
%       NA - Numerical aperature of the element
%
%       Manufacture - The element manufacturer
%
%       ElementID - Element ID from the manufacturer
%
%       Price - Price of the element
%
%       DateOfInformation - The date this information is from. Empty if the
%       date is unknown. Class string.
%
% opticalElement3D Methods:
%   
%       drawElement - Create the graphical objects necessary to plot
%       the element
%
%       transformRays - Transform rays from the global coordinate system to
%       the element's coordinate system
%
%       invTransformRays - Transform rays from the element's coordinate
%       system to the global coordinate system
%
% See also OPTICALELEMENT_SURFACE SPHERICALSURFACE ASPHERICALSURFACE
% PLANESURFACE CYLINDRICALSURFACE CONESURFACE APERTURESURFACE
% PARABOLICSURFACE BACKSURFACE

% James Kapaldo, January 2016 

% 2016, 01, 23 -- Added in sequentialRayTrace property, ExtentAA property,
% Matrix property, removeRaysOutsideCrossSection method, transformRays
% method, and invTransformRays method. Modified updatePosition method to
% just use the Orientation and Origin for transformation, instead of trying
% to calculate vectors to rotate about and angles to rotate by. Modified
% drawCube so that min vertex is column 1 and max vertex is column 8. Added
% in function createAxis and am not drawing an axis on the optical elements
% with the 'Tag', 'cordSys' (default is visible = off). Added in set
% methods from crossSectionSize and sequentialRayTrace. Modfied the rayDat
% array that is being used: removed theta, phi, and I (intensity); rayDat
% is now nx8 with the 8 columns being x,y,z,uv_x,uv_y,uv_z,w,l. 

    properties(SetAccess = public, GetAccess = public)        
        Type                = []; % Type of the optical element: {lens, spacial filter, mirror}
        SubType             = [];
        Focal               = []; % focal lengths, normally at 1x2 vector with (focal, back focal)
        NA                  = []; % Numerical aperature
    end
    
    properties(SetAccess = public, GetAccess = public)
        crossSection        = [];
        crossSectionSize    = [];           % If crossSection = 'circle', then numel(crossSectionSize) = 1 and the value should be the diameter. If crossSection = 'rectangle', then, if numel = 1, the shape will be a square with edge length Dimension; if numel = 2, the shape is a rectangle with an edge length crossSectionSize(1) along Orientation(:,2) and an edge length crossSectionSize(2) along Orienation(:,3)
        Orientation         = eye(3);       % size = [3,3], Orientation(:,1) is a unit vector pointing along the optical axis (center axis of a lens). Orientation(:,2) is a vector perpendicular to Orientation(:,1). Orientation(:,3) is a unit vector along cross(Orientation(:,1),Orientation(:,2)). -- When creating the element, only Orientation(:,1) and Orientation(:,2) are strictly required. If Orientation(:,3) is not given, then it will be calculated from Orientation(:,1:2).
        Origin              = zeros(3,1);   % size = [3,1], vector giving the center of the element.
        Surfaces            = nullSurface;  % heterogeneous array of class opticalElement_Surface with size [1 n]'.
        
        Material            = [];           % cell array of strings giving the materials in the order they appear from the Surfaces array
    end
    
    properties(SetAccess = protected)
        Extent           = []; % This is a derived array of 8 points giving the bounding volume of the element
        ExtentAA         = []; % The axis alinged extent of the element. size = [2,3], first row is minimum point, second row is maximum point.
        Matrix           = [];
    end
    
    properties(SetAccess = public, GetAccess = public)
        ShapeGroup       = matlab.graphics.primitive.Transform;
%         rayTraceData     = [];
    end
    
    properties(SetAccess = public, GetAccess = public)        
        Manufacture         = [];
        ElementID           = [];
        Price               = [];
        DateOfInformation   = [];
    end
    
    properties (Hidden = true)
        sequentialRayTrace  = true;         % set if ray tracing is done sequentially through the surfaces, or non-sequentially.
        ShapeGroupPrimary = matlab.graphics.primitive.Transform;
        ExtentPrimary = [];
    end
    
%     
%     properties(SetAccess = {?opticalElement, ?opticalSystem}, ...
%                GetAccess = public, Hidden = false)
%            
%         wasModifiedListener = [];
%         
%         wasModified_Ray     = []; %boolian field set to true if the element was moved or modified in some why after the last ray trace or scalar diffraction calculation.   
%         
%         Next            = []; %next element
%         Prev            = []; %previous element
%         Parent          = []; %opticalSystem that the element belongs to
%         MovementBounds  = []; %Range which the element may be moved.
%         
%         plottedRays         = [];
%     end

    
    methods 
        
%         function el = opticalElement(varargin)
%             
%         end %OpticalElement

        function set.Orientation(obj,value)
            validateattributes(value,{'double'},{'finite','nrows',3})
            if size(value,2) < 2 || size(value,2) > 3
                error('opticalElement3D:OrientationWrongSize','The size of Orientation must be at least 3x2 and at most 3x3')
            end
            
            if size(value,2) == 2
                if ~(all(abs(sum(value.^2) - 1) < 10*eps) && abs(sum(value(:,1).*value(:,2))) < 10*eps)
                    error('opticalElement3D:OrientationWrongSize','The columns of Orientation must be unit vectors and the columns must be perpendicular to each other')
                end
                value(:,end+1) = cross(value(:,1),value(:,2));
                value(:,end) = value(:,end)/sqrt(sum(value(:,end).^2));
            else % 3 columns were given
                if det(value) ~= 1 
                    error('opticalElement3D:OrientationWrongSize','The columns of Orientation must be unit vectors and the columns must be perpendicular to each other and form a right handed coordinate system.')
                end
            end
            
            obj.Orientation = value;
            if ~isempty(obj.ExtentPrimary) %#ok<MCSUP>
                updatePosition(obj);
            end
        end
        
        function set.Origin(obj, value)
            validateattributes(value,{'double'},{'finite','numel',3})
            obj.Origin = value(:);
            if ~isempty(obj.ExtentPrimary) %#ok<MCSUP>
                updatePosition(obj);
            end
        end

        function set.sequentialRayTrace(obj, value)
            validateattributes(value,{'logical'},{'numel',1})
            obj.sequentialRayTrace = value;
        end
        
        function set.crossSectionSize(obj, value)
            if numel(value)>2 || numel(value)<1 || ~all(isfinite(value))
                error('setCrossSectionSize:badInput','While setting crossSectionSize, expected input to have 1 or 2 finite elements.');
            end
            
            obj.crossSectionSize = value(:)';
        end
        
        function set.crossSection(obj, value)
            value = validatestring(value,{'circle','rectangle'});
            obj.crossSection = value;
        end
        
        function delete(el)
            % delete  Delete the optical element.
            
            % First need to remove the element from any axis if it has been
            % plotted.
            el.ShapeGroup.Parent = matlab.graphics.GraphicsPlaceholder.empty;
        end

        function resetDefaults(obj)
            % resetDefaults  clear the optical element and set all
            % properties to the default values.
            if isvalid(o)
                mc = metaclass(obj);
                mp = mc.PropertyList;
                for k=1:length(mp)
                    if mp(k).HasDefault && ~strcmp(mp(k).SetAccess,'private')
                        obj.(mp(k).Name) = mp(k).DefaultValue;
                    end
                end
            end
        end

        
%         function empty = isempty(el)
%             
%             if ~isvalid(el)
%                 fprintf(2, 'Deleted object.\n\n')
%                 empty = 1;
%                 return;
%             end
%             
%             if isscalar(el)
%                 if isempty(el.Type)
%                     empty = 1;
%                 else
%                     empty = 0;
%                 end
%             else
%                 for i = 1:numel(el)
%                     if ~isempty(el(i).Type)
%                         empty = 0;
%                         return;
%                     end
%                 end
%                 empty = 1;
%             end
%         end
        
%         function disp(el, varargin)
%             % varargin may hold an optional index parameter that marks the
%             % position of el inside an opticalSystem
%             
%             if ~isvalid(el)
%                 fprintf(2, 'Deleted object.\n\n') %#ok<*PRTCAL>
%                 return;
%             end
%             
%             if (isscalar(el))
%                 if isempty(el.Type)
%                     fprintf('Empty optical element.\n\n')
%                 else
%                     if nargin > 1
%                         fprintf('Optical element in optical system, %s, at position %d: \n', ...
%                             varargin{1}, varargin{2})
%                     else
%                         fprintf('Optical element: \n')
%                     end
%                     fprintf('\t Type: \t\t\t\t%-29s\n', el.Type)
%                     if strcmp(el.Type, 'Axicon')
%                         param = -el.SimpleSurface{2}{2};
%                         str = 'Angle: \t';
%                         unit = ['[' char(176) ']'];
%                     elseif strcmp(el.Type, 'Spacial Filter')
%                         param = NaN;
%                     else
%                         param = el.Focal(1);
%                         str = 'Focal length: ';
%                         unit = '[mm]';
%                     end
%                     
%                     if ~isnan(param) %print out this line of information for all elements but spacial filters
%                         fprintf(['\t ' str '\t\t%- 16.2f %s\n'], param, unit)   
%                     end
%                     
%                     fprintf('\t Width: \t\t\t%- 16.2f [mm]\n', el.Width)
%                     fprintf('\t Z: \t\t\t\t%- 16.2f [mm]\n', el.Z)
%                     if ~isempty(el.MovementBounds)
%                         fprintf('\t Movement Bounds: \t[%-6.2f, %-6.2f] [mm]\n\n', ...
%                             el.MovementBounds(1), el.MovementBounds(2))
%                     end
%                 end
%             else
%                 dims = size(el);
%                 dims = num2str(dims);
%                 dims(3:3:length(dims)) = [];
%                 dims(2:2:length(dims)) = 'x';
%                 fprintf('%s array of optical elements\n\n', dims)
%             end
%         end %disp
        
%         function get(el)
%             builtin('disp',el)
%         end
        
%         flipElement(el);
        
    end %methods
    
    methods(Access = public)%{?opticalElement, ?opticalSystem}, Hidden = true)
        
        function drawElement(el)
            % drawElement  Create the graphical objects necessary to plot
            % the element
            %
            % This function will create a transform object and save it in
            % the ShapeGroup property of the element. The Parent of this
            % transform group will be set to empty. To display the element
            % in an axis, set the Parent to the axes. 
            %   e.g.: el.ShapeGroup.Parent = gca
            
            % This function will actually just copy the already created
            % transform group from ShapeGroupPrimary to ShapeGroup. This is
            % done so that when a figure is closed that held the
            % ShapeGroup, it is not necessary to recreate all of the
            % element surfaces. We just create the element group once and
            % store it in ShapeGroupPrimary, and then we copy to ShapeGroup
            % when drawElement is called.
            try
                el.ShapeGroup = copy(el.ShapeGroupPrimary);
                el.ShapeGroup.DeleteFcn = {@shapeGroupBeingDeleted, el};
            catch

            end
            function shapeGroupBeingDeleted(~,~, el)
                try
                    el.ShapeGroup = copy(el.ShapeGroupPrimary);
                catch
                end
            end
        end
        
        function rayDat = transformRays(el,rayDat)
            % transformRays  Transform rays from the global coordinate
            % system to the element's coordinate system.
            %
            % Inputs:
            %   el - The optical element
            %
            %   rayDat - An array of ray data (n x 8) where n is the number
            %   of rays. Columns 1:3 must be the ray's positions, and
            %   columns 4:6 must be the ray's directions.
            %
            % Outputs:
            %   rayDat - The transformed ray data
            %
            % See also INVTRANSFORMRAYS
            
            one = ones(1,size(rayDat,1));
            
            if size(rayDat,3) == 1
                temp = el.Matrix \ [rayDat(:,1:3)'; one];
                rayDat(:,1:3) = temp(1:3,:)';

                rayDat(:,4:6) = (el.Orientation' * rayDat(:,4:6)')';
            else
                for i = 1:size(rayDat,3)
                    temp = el.Matrix \ [rayDat(:,1:3,i)'; one];
                    rayDat(:,1:3,i) = temp(1:3,:)';

                    rayDat(:,4:6,i) = (el.Orientation' * rayDat(:,4:6,i)')';
                end
            end
        end
        
        function rayDat = invTransformRays(el,rayDat)
            % invTransformRays  Transform rays from the element's
            % coordinate system to the global coordinate system.
            %
            % Inputs:
            %   el - The optical element
            %
            %   rayDat - An array of ray data (n x 8) where n is the number
            %   of rays. Columns 1:3 must be the ray's positions, and
            %   columns 4:6 must be the ray's directions.
            %
            % Outputs:
            %   rayDat - The transformed ray data
            %
            % See also TRANSFORMRAYS
            
            one = ones(1,size(rayDat,1));
            
            if size(rayDat,3) == 1
                temp = el.Matrix * [rayDat(:,1:3)'; one];
                rayDat(:,1:3) = temp(1:3,:)';

                rayDat(:,4:6) = (el.Orientation * rayDat(:,4:6)')';
            else
                for i = 1:size(rayDat,3)
                    temp = el.Matrix * [rayDat(:,1:3,i)'; one];
                    rayDat(:,1:3,i) = temp(1:3,:)';

                    rayDat(:,4:6,i) = (el.Orientation * rayDat(:,4:6,i)')';
                end
            end
        end
        
%         createLens(el, varargin);
%         createSpacialFilter(el, varargin);
%         
%         [rIn, data, phaseDist] = rayTraceSimpleElement(el, rIn);

    end
    
    methods (Hidden = true) % Access = protected
        
        function updatePosition(el)
            % updatePosition  Update the location and orientation of the
            % element using the Orientation and Origin properties.
            
            % This function will be called after Orientation or Center has
            % been updated.
            
            rot = [el.Orientation, el.Origin; zeros(1,3), 1];
            el.ShapeGroupPrimary.Matrix = rot;
            el.Matrix = rot;
            try
                el.ShapeGroup.Matrix = rot;
            catch
            end
            
            
            v = rot*[el.ExtentPrimary.Vertices;ones(1,8)];
            n = rot*[el.ExtentPrimary.Normals;ones(1,6)];
            
            el.Extent.Vertices = v(1:3,:);
            el.Extent.Normals = n(1:3,:);
            
            el.ExtentAA = [min(el.Extent.Vertices,[],2), max(el.Extent.Vertices,[],2)]';
        end
        
        function drawElementFirstTime(el)
            % drawElement Construct a 3D model of an optical element(s)
            %   This function will take in an array of type
            %   opticalElement3D and create a 3D model for each of them.
            %   The models will be Transform groups stored in the property
            %   ShapeGroup. Be default the Transform groups do not have a
            %   Parent; to display the models, set the Parent property of
            %   each ShapeGroup.
            
            for j = 1:numel(el)
                el(j).ShapeGroupPrimary = matlab.graphics.primitive.Transform;               
                
                % First, create all but the aperatureSurfaces, since
                % aperatureSurfaces get in the way of the Material indexing
                realSurfs = find(~isa(el(j).Surfaces,'apertureSurface'));
                
                createSurfPatch(el(j).Surfaces(realSurfs(1)),el(j),realSurfs(1),el(j).Material{1})

                for i = 2:numel(realSurfs)
                    createSurfPatch(el(j).Surfaces(realSurfs(i)),el(j),realSurfs(i),el(j).Material{i-1})

                    edgeData_x = [el(j).Surfaces(realSurfs(i-1)).boundingEdge(:,1)';el(j).Surfaces(realSurfs(i)).boundingEdge(:,1)'];
                    edgeData_y = [el(j).Surfaces(realSurfs(i-1)).boundingEdge(:,2)';el(j).Surfaces(realSurfs(i)).boundingEdge(:,2)'];
                    edgeData_z = [el(j).Surfaces(realSurfs(i-1)).boundingEdge(:,3)';el(j).Surfaces(realSurfs(i)).boundingEdge(:,3)'];

%                     edgeLength = abs(edgeData_x(1,:) - edgeData_x(2,:));
%                     fprintf('Min, max edge length for edge %d: [%0.2f %0.2f]\n',i-1,min(edgeLength),max(edgeLength));
                    
                    [f,v] = surf2patch(edgeData_x,edgeData_y,edgeData_z);

                    [v,f] = opticalElement_Surface.removeDubs(v,f);
                    
                    options = opticalElement_Surface.displayOptions(el(j).Surfaces(realSurfs(i)), el(j).Material{i-1});
            
                    patch('Vertices', v, ...
                          'Faces', f, ...
                          'Parent', el(j).ShapeGroupPrimary, ...
                          options);
                end

                % Now, creat the aperature surfaces
                aprtrSurfs = find(isa(el(j).Surfaces,'apertureSurface'));
                
                for i = 1:numel(aprtrSurfs)
                    createSurfPatch(el(j).Surfaces(aprtrSurfs(i)),el(j),aprtrSurfs(i),'none')
                end
                
                % Plot an axis
                h = createAxis(eye(3),zeros(3,1),5);
                set(h,'Parent',el(j).ShapeGroupPrimary);
                
                % Get the bounding box for the element, the extent
                extent = [inf(3,1),-inf(3,1)];
                
                for i = 1:numel(el(j).ShapeGroupPrimary.Children)
                    child = el(j).ShapeGroupPrimary.Children(i);
                    if strcmp(child.Type,'hggroup')
                        for k = 1:numel(child.Children)
                            if strcmp(child.Children(k).Type,'patch')
                                verts = child.Children(k).Vertices;
                                if size(verts,2) == 3
                                    verts = verts';
                                end

                                extent(:,1) = min(extent(:,1),min(verts,[],2));
                                extent(:,2) = max(extent(:,2),max(verts,[],2));
                            end
                        end
                    elseif strcmp(child.Type,'patch')
                        verts = child.Vertices;
                        if size(verts,2) == 3
                            verts = verts';
                        end

                        extent(:,1) = min(extent(:,1),min(verts,[],2));
                        extent(:,2) = max(extent(:,2),max(verts,[],2));
                    end
                end
                
                center = sum(extent,2)/2;

                edgeLengths = abs(extent(:,1)-extent(:,2));
                [v,f,n] = drawCube(center,edgeLengths);

                patch('Vertices', v', ...
                      'Faces', f, ...
                      'Parent', el(j).ShapeGroupPrimary, ...
                      'FaceColor','none',...
                      'EdgeColor','k',...
                      'LineStyle','--',...
                      'Tag','BoundingBox');
                  
                el(j).ExtentPrimary.Vertices = v;
                el(j).ExtentPrimary.Faces = f;
                el(j).ExtentPrimary.Normals = n;
                
                updatePosition(el(j));
                
            end
            
            drawElement(el)
        
        end
        
        function [rayDat, inside] = removeRaysOutsideCrossSection(el,rayDat)
            % removeRaysOutsideCrossSection  (pretty self explanitory what
            % the function does)
            %
            % This function should be called after calling goToSurf on a
            % surface to remove points that fall outside of the elements
            % crossSection. 
            %
            % Inputs:
            %   el - opticalElement3D
            %   rayDat - ray data
            %
            % Outputs:
            %   rayDat - the ray data witth positions that fall outside of
            %   the cross section replaced with NaN's
            %
            %   inside - A logical array of size [size(rayDat,1),1] whose
            %   elements are true if the ray is inside of the crossSection
            %   and false otherwise.
            
            
            if size(size(rayDat)) ~= 2
                error('removeRaysOutsideCrossSection:badRayDataSize','The ray data may not have a 3rd dimension for this function.')
            end
            
            switch el.crossSection
                case 'circle'
                    inside = sqrt(sum(rayDat(:,2:3).^2,2)) < el.crossSectionSize(1);                     
                case 'rectangle'
                    if numel(el.crossSectionSize) == 1
                        inside = abs(rayDat(:,2)) < el.crossSectionSizes/2 & abs(rayDat(:,3)) < el.crossSectionSizes/2;
                    else
                        inside = abs(rayDat(:,2)) < el.crossSecitonSizes(1)/2 & abs(rayDat(:,3)) < el.crossSecitonSizes(2)/2;
                    end
                otherwise
                    error('removeRaysOutsideCrossSection:badCrossSection','An unknown crossSection shape, %s, has occured.', el.crossSection)
            end
            
            if size(rayDat,2) >= 8 % do this because sometimes I will only call this function with the ray positions rayDat(:,1:3) and not with the entire rayDat array
                rayDat(~inside,[1:3,8]) = NaN;
            else
                rayDat(~inside,1:3) = NaN;
            end
            
        end
    end
    
%     methods(Access = public, ... %{?opticalElement, ?opticalSystem}, Hidden = true, ...
%             Static = true)
%         
%         rIn = rayTraceFreeSpace(rIn,d,n);
%         rIn = conicalSurface(rIn, A, n1, n2);
%         rIn = sphericalSurface(rIn, R, n1, n2);
%         rIn = asphericalSurface(rIn, params, n1, n2);
%         
%         h = distanceToConicalSurf(rIn, A, y, n);
%         h = distanceToSphericalSurf(rIn, R, y, n);
%         h = distanceToAsphericalSurf(rIn, params, y, n, range);
%         n = refractiveIndex(material, x);
%         
%         simpleSurf = constructSimpleSurfaceArray(names, params, varargin);
%     end
    
end %classdef

function [v, f, n] = drawCube ( origin, size )
% drawCube  Create vertices and faces for rectangular cuboid with size
% "size" and origin "origin".
    v = [0 0 0;
         0 1 0;
         1 1 0;
         1 0 0;
         0 0 1;
         0 1 1;
         1 0 1;
         1 1 1]'-0.5;
    
    
    
    for i = 1:3
        v(i,:) = v(i,:)*size(i) + origin(i);
    end
    f = [3 4 7 8;2 3 8 6;5 6 8 7;1 2 6 5;5 7 4 1;4 3 2 1];
    
    n = zeros(3,6);
    for i = 1:6
        n(:,i) = cross(v(:,f(i,2))-v(:,f(i,3)), v(:,f(i,1))-v(:,f(i,2)));
        n(:,i) = n(:,i)/norm(n(:,i));
    end
    
end

function h = createAxis(X, O, length_v)

    c = [1      0.3137 0.2784; %red
         1      0.8431 0;      %yellow
         0.1804 0.5451 0.3412];%green

for i = 3:-1:1
    x = O(1)+[0, X(1,i)]*length_v;
    y = O(2)+[0, X(2,i)]*length_v;
    z = O(3)+[0, X(3,i)]*length_v;
    h(i) = matlab.graphics.primitive.Line;
    h(i).XData = x;
    h(i).YData = y;
    h(i).ZData = z;
    h(i).Color = c(i,:);
    h(i).LineWidth = 2;
    h(i).Tag = 'cordSys';
    h(i).Visible = 'off';
end
end
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
