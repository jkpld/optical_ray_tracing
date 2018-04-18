%% Example ray tracing
% Ray trace a point light source being colminated by a parabolic mirror and
% then focues by an aspheric lens.

%% Load optics
try delete(optic(1)), catch, end
try delete(optic(2)), catch, end
optic(2) = readLensFile_json('lens_files\RC12SMA_F01.json');
optic(1) = readLensFile_json('lens_files\eo48536.json');

%% Position optics

% Position the aspheric lens so that it focuses light to the origin.
opAx = [-1;0;0];
opAx = opAx/norm(opAx);
uv2 = gramSchmidt1([0;1;0]',opAx')';
uv2 = uv2/norm(uv2);

optic(1).Orientation = [opAx, uv2];
optic(1).Origin = [-optic(1).Focal(2) - range(optic(1).ExtentPrimary.Vertices(1,:)), 0, 0];

% Position the parabolic
optic(2).Orientation = [1 0; 0 0; 0 -1];
optic(2).Origin = [-140;0;0];

% Plot the optical system
fig = figure('pos',[496, 394, 823, 363]);
ax = axes;

for i=1:numel(optic)
    try
        optic(i).ShapeGroup.Parent = ax; 
    catch
        % If the optics were already placed in a figure, and that figure
        % was closed, then the ShapeGroup would have been destroyed.
        % Therefore, recreat the ShapeGroup's for all of the optics.
        drawElementFirstTime(optic);
        optic(i).ShapeGroup.Parent = ax; 
    end
end

daspect([1 1 1]);
setTheme(fig,'light')
ax.XLim = [-170, 20];
ax.XLabel.String = "Optical axis (mm)";
view(26,13)

%% Initialize rays

% Build point light source directions
n = icoSphere(4); % create normal vectors for point source distributed around a sphere
n(n(:,3)<0,:) = []; % remove all rays pointing in the positive x direction

% Initialize rays
rays = zeros(size(n,1), 8,1);
rays(:,1) = -140; % x pos
rays(:,2) = 0; % y pos
rays(:,3) = -50.8; % z pos
rays(:,4:6) = n; % ray direction
rays(:,7) = 587.6; % ray wavelength

% Remove rays that do not intersect first optic
intersects = rayIntersectAABB(rays, optic(2).ExtentAA);
rays(~intersects,:) = [];

%% Ray trace elements
% Ray trace both elements and then propogate rays an extra 30.
raysO1 = rayTraceElement(optic(2), rays);
raysO2 = rayTraceElement(optic(1), raysO1(:,:,end));
raysO2(:,:,end+1) =  propagation(raysO2(:,:,end),30,1);

% Combine all rays and plot
rayDat = cat(3, rays, raysO1, raysO2);

% Plot results
try delete(lines234), catch, end
lines234 = line(squeeze(rayDat(:,1,:))', squeeze(rayDat(:,2,:))', squeeze(rayDat(:,3,:))','Color','k','Marker','.','MarkerEdgeColor','r','Parent',ax);