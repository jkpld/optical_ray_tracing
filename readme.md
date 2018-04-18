# Ray tracing
Methods for performing 3d sequential ray tracing of optical systems.

Each optical element is made up from a set of optical surfaces, and each surface contains the methods for determining the reflection/refraction of rays off of the surface.

![optical_system_w_rays](/examples/docs/optical_system_w_rays.png)

## Example usage 1
Ray trace a point light source being colminated by a parabolic mirror and then focues by an aspheric lens.

### Load optics
```Matlab
try delete(optic(1)), catch, end
try delete(optic(2)), catch, end
optic(2) = readLensFile_json('lens_files\RC12SMA_F01.json');
optic(1) = readLensFile_json('lens_files\eo48536.json');
```

### Position optics
Position the aspheric lens so that it focuses light to the origin.
```Matlab
opAx = [-1;0;0];
opAx = opAx/norm(opAx);
uv2 = gramSchmidt1([0;1;0]',opAx')';
uv2 = uv2/norm(uv2);

optic(1).Orientation = [opAx, uv2];
optic(1).Origin = [-optic(1).Focal(2) - range(optic(1).ExtentPrimary.Vertices(1,:)), 0, 0];
```

Position the parabolic
```Matlab
optic(2).Orientation = [1 0; 0 0; 0 -1];
optic(2).Origin = [-140;0;0];
```

Plot the optical system
```Matlab
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
```
![optical_system](/examples/docs/optical_system.png)

### Initialize rays
Build point light source directions
```Matlab
n = icoSphere(4); % create normal vectors for point source distributed around a sphere
n(n(:,3)<0,:) = []; % remove all rays pointing in the positive x direction
```

Initialize rays
```Matlab
rays = zeros(size(n,1), 8,1);
rays(:,1) = -140; % x pos
rays(:,2) = 0; % y pos
rays(:,3) = -50.8; % z pos
rays(:,4:6) = n; % ray direction
rays(:,7) = 587.6; % ray wavelength
```

Remove rays that do not intersect first optic
```Matlab
intersects = rayIntersectAABB(rays, optic(2).ExtentAA);
rays(~intersects,:) = [];
```

### Ray trace elements
Ray trace both elements and then propagate rays an extra 30.
```Matlab
raysO1 = rayTraceElement(optic(2), rays);
raysO2 = rayTraceElement(optic(1), raysO1(:,:,end));
raysO2(:,:,end+1) =  propagation(raysO2(:,:,end),30,1);
```

Combine all rays and plot
```Matlab
rayDat = cat(3, rays, raysO1, raysO2);
```

Plot results
```Matlab
try delete(lines234), catch, end
lines234 = line(squeeze(rayDat(:,1,:))', squeeze(rayDat(:,2,:))', squeeze(rayDat(:,3,:))','Color','k','Marker','.','MarkerEdgeColor','r','Parent',ax);
```
![optical_system_w_rays](/examples/docs/optical_system_w_rays.png)


## Example Usage 2
Plot the focal length shift of an aspheric lens, made with fused silica, as the wavelength changes

### Create optical element
```Matlab
try delete(optic(1)), catch, end
optic(1) = readLensFile_json('K:\sequential_ray_tracing\lens_files\eo48536.json');
```

Orient element
```Matlab
opAx = [-1;0.0;0.0];
opAx = opAx/norm(opAx);
uv2 = gramSchmidt1([0;1;0]',opAx')';
uv2 = uv2/norm(uv2);

optic(1).Orientation = [opAx, uv2];
optic(1).Origin = [-optic(1).Focal(2) - range(optic(1).ExtentPrimary.Vertices(1,:)), 0, 0];
```

### Create initial ray data
Create a ring of rays each with a different wavelength going from 210 to 880 nm.
```matlab
N = 200; % number of rays
r = 7; % radius of rays
th = linspace(0,2*pi,N).'; % azimuthal angle
r0 = [-40*ones(N,1),cos(th)*r,sin(th)*r];
uv0 = zeros(N,3); % directions
uv0(:,1) = 1;
w = linspace(210,880,N).'; % wavelengths
rays = [r0,uv0,w,zeros(N,1)];
```

Remove any rays that do not intersect the first object.
```matlab
intersects = rayIntersectAABB(rays, optic(1).ExtentAA);
rays(~intersects,:) = [];
```

Propagate the rays through the lens and to the optical axis
```matlab
rayOut = rayTraceElement(optic(1), rays);
```

Propagate rays to optical axis
```matlab
d = -rayOut(:,2,end)./rayOut(:,5,end);
rayDat = cat(3,rays,rayOut,propagation(rayOut(:,:,end),d,1));
```

Create figure with lens and rays
```matlab
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
```
![focal_shift1](/examples/docs/focal_shift_1.png)

### Plot focal shift
```matlab
figure
line(rayDat(:,1,end), w,'color','r','linewidth',2);
axis tight
xlabel('Focal shift / mm')
ylabel('Wavelength / nm')
title(sprintf('Focal shift of asphere (EO: %s)',optic.ElementID))
setTheme(gcf,'light')
ylim([200,900])
```
![focal_shift2](/examples/docs/focal_shift_2.png)

## Code structure
Each optical element is stored in a class `opticalElement3D` that contains information about the element, its location and orientation in 3D space, and the set of surfaces and materials that make up the element. Each of the surfaces is a concrete implementation of the abstract class `opticalElement_Surface`. These surface classes contain methods for propogating rays to the surface and reflecting or refracting off of the surface. The rays polarization can also be considered.
