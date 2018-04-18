%% Focal shift of EO48536
% Determine the focal shift of an aspheric lens

%% Create optical element
try delete(optic(1)), catch, end
optic(1) = readLensFile_json('K:\sequential_ray_tracing\lens_files\eo48536.json');

% Orient element
opAx = [-1;0.0;0.0];
opAx = opAx/norm(opAx);
uv2 = gramSchmidt1([0;1;0]',opAx')';
uv2 = uv2/norm(uv2);

optic(1).Orientation = [opAx, uv2];
optic(1).Origin = [-optic(1).Focal(2) - range(optic(1).ExtentPrimary.Vertices(1,:)), 0, 0];

%% Create initial ray data - ring of rays each with a different wavelength
N = 200;
r = 7;
th = linspace(0,2*pi,N).';
r0 = [-40*ones(N,1),cos(th)*r,sin(th)*r];
uv0 = zeros(N,3);
uv0(:,1) = 1;
w = linspace(210,880,N).';
rays = [r0,uv0,w,zeros(N,1)];

% Remove any rays that do not intersect the first object.
intersects = rayIntersectAABB(rays, optic(1).ExtentAA);
rays(~intersects,:) = [];

% Propogate the rays through the lens and to the optical axis
rayOut = rayTraceElement(optic(1), rays);

% Propogate rays to optical axis
d = -rayOut(:,2,end)./rayOut(:,5,end);
rayDat = cat(3,rays,rayOut,propagation(rayOut(:,:,end),d,1));

% Create figure with lens and rays
figure
axLens = axes;
optic(1).ShapeGroup.Parent = axLens;
rys = @(i) [squeeze(rayDat(:,i,:))'; nan(1,N)];
idx = @(x,varargin) x(varargin{:});
fl = @(x) x(:);
pltRys = @(i) fl(idx(rys(i),':',1:7:N));

lines234 = line(pltRys(1),pltRys(2),pltRys(3),'Color','k','Marker','.','MarkerEdgeColor','r','Parent',axLens,'LineStyle','-','LineWidth',0.5);
lines234.Color(4) = 0.2;

xlabel('Optical axis / mm')
ylabel('y / mm')
zlabel('z / mm')
setTheme(gcf,'light')
daspect([1 1 1]);
axLens.Visible = 'off';
axis tight
optic(1).ShapeGroup.Visible = 'on';
view(0,0)


%% Plot focal shift
figure
line(rayDat(:,1,end), w,'color','r','linewidth',2);
axis tight
xlabel('Focal shift / mm')
ylabel('Wavelength / nm')
title(sprintf('Focal shift of asphere (EO: %s)',optic.ElementID))
setTheme(gcf,'light')
ylim([200,900])