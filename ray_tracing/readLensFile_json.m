function el = readLensFile_json(file)
el = opticalElement3D;
el.Focal = zeros(1,2);

txt = fileread(file);
dat2 = jsondecode(txt);

fn = fieldnames(dat2);
for i = 1:numel(fn)
    dat.(lower(fn{i})) = dat2.(fn{i});
    if strcmpi(fn{i},'surfaces')
        for j = 1:numel(dat2.(fn{i}))
            fn2 = fieldnames(dat2.(fn{i}){j});
            dat3 = struct();

            for k = 1:numel(fn2)
                dat3.(lower(fn2{k})) = dat2.(fn{i}){j}.(fn2{k});
            end
            dat.(lower(fn{i})){j} = dat3;
        end
    end
end

fn = string(fieldnames(dat));

for i = fn'
    switch i
        case lower("Manufacturer")
            el.Manufacture = dat.(char(i));
        case lower("ID")
            el.ElementID = dat.(char(i));
        case lower("Type")
            el.Type = lower(dat.(char(i)));
        case lower("SubType")
            el.SubType = lower(dat.(char(i)));
        case lower("crossSection")
            if ~any(string(lower(dat.(char(i)))) == ["circle", "rectangle"])
                error('readLensFile:badCrossSection','The crossSection must be either "circle" or "rectangle".');
            end
            el.crossSection = lower(dat.(char(i)));
        case lower("crossSectionSize")
            temp = dat.(char(i));
            if numel(temp)>2 || numel(temp)<1 || ~all(isfinite(temp))
                error('readLensFile:badCrossSectionSize','The crossSectionSize must have 1 or 2 finite elements.');
            end

            if strcmp(el.crossSection,'circle')
                temp = temp/2;
            end

            el.crossSectionSize = temp(:)';
        case lower('sequential')
            temp = dat.(char(i));

            if ~(numel(temp)==1 && (temp ==0 || temp ==1))
                error('readLensFile:badSequential','Bad sequential, must be 0 or 1');
            end

            el.sequentialRayTrace = logical(temp);
        case lower("clearAperture")
            temp =  dat.(char(i));

            if isempty(temp)
                clearAperture = [];
            else
                if strcmp(el.crossSection,'circle')
                    temp = temp/2;
                end
                if numel(temp) ~= numel(el.crossSectionSize) || ~all(isfinite(temp)) || ~all(temp<=el.crossSectionSize)
                    error('readLensFile:badclearAperture','Bad clearAperture, number of elements must be the same size as crossSectionSize, the values must be finite and less than the crossSectionSize');
                end
                clearAperture = temp(:)';
            end
        case lower("f")
            temp =  dat.(char(i));
            if ~isempty(temp)
                if numel(temp) ~= 1 || ~isfinite(temp)
                    error('readLensFile:badFocalLength','Bad focal length');
                end
            end
            el.Focal(1) = temp;
        case lower("fb")
            temp = dat.(char(i));
            if ~isempty(temp)
                if numel(temp) ~= 1 || ~isfinite(temp)
                    error('readLensFile:badbackFocalLength','Bad back focal length');
                end
                el.Focal(2) = temp;
            else
                el.Focal(2) = [];
            end
        case lower("NA")
            temp = dat.(char(i));
            if ~isempty(temp)
                if numel(temp) ~= 1 || ~isfinite(temp)
                    error('readLensFile:badNA','Bad NA value');
                end
            end
            el.NA = temp;
        case lower("Price")
            temp = dat.(char(i));
            el.Price = temp(1);
        case lower("DateOfInformation")
            temp = dat.(char(i));
            el.DateOfInformation = temp;
        case lower("edgeThickness")
            temp = dat.(char(i));
            %             if numel(temp) ~= 1 %|| ~isfinite(temp)
            %                 error('readLensFile:badedgeThickness','Bad edgeThickness');
            %             end
            edgeThickness = temp;
        case lower("Surfaces")
            numSurfaces = numel(dat.surfaces);
            for j = 1:numSurfaces
                fn_s = string(fieldnames(dat.surfaces{j}));

                for k = fn_s'
                    switch k
                        case lower("material")
                            if j ~= numSurfaces
                                el.Material{length(el.Material)+1} = dat.surfaces{j}.(char(k));
                            end
                        case lower("surfaceType")
                            surfType = dat.surfaces{j}.(char(k));
                            switch surfType
                                case 'sphere'
                                    el.Surfaces(j) = sphericalSurface;
                                    el.Surfaces(j).clearAperture = clearAperture;

                                    el = assignSurfProp(el,'R',dat,'R', j);
                                    el = assignSurfProp(el,'Center',dat,'apex', j);

                                    el.Surfaces(j).Center = el.Surfaces(j).Center - [dat.surfaces{j}.R; 0; 0];
                                case 'plane'
                                    el.Surfaces(j) = planeSurface;
                                    el.Surfaces(j).clearAperture = clearAperture;

                                    el = assignSurfProp(el,'normVec',dat,'normVec', j);
                                    el = assignSurfProp(el,'Point',dat,'point', j);

                                case 'backSurface'
                                    el.Surfaces(j) = backSurface;
                                    el.Surfaces(j).clearAperture = [];

                                    el = assignSurfProp(el,'normVec',dat,'normVec', j);
                                    el = assignSurfProp(el,'Point',dat,'point', j);

                                case 'parabola'
                                    el.Surfaces(j) = parabolicSurface;
                                    el.Surfaces(j).clearAperture = clearAperture;

                                    el = assignSurfProp(el,'PFL',dat,'PFL', j);
                                    el = assignSurfProp(el,'EFL',dat,'EFL', j);
                                    el = assignSurfProp(el,'Angle',dat,'Angle', j);
                                    el = assignSurfProp(el,'baseDist',dat,'baseDist', j);

                                case 'aperture'

                                    el.Surfaces(j) = apertureSurface;
                                    el.Surfaces(j).clearAperture = clearAperture;

                                    el = assignSurfProp(el,'Shape',dat,'shape', j);
                                    el = assignSurfProp(el,'normVec',dat,'normVec', j);
                                    el = assignSurfProp(el,'Center',dat,'center', j);
                                    el = assignSurfProp(el,'Extent',dat,'Extent', j);

                                    el.Material(end) = [];

                                case 'asphere'
                                    el.Surfaces(j) = asphericalSurface;

                                    el = assignSurfProp(el,'Center',dat,'apex', j);
                                    el = assignSurfProp(el,'R',dat,'R', j);
                                    el = assignSurfProp(el,'k',dat,'conic', j);
                                    el = assignSurfProp(el,'A',dat,'coeffs', j);

                                otherwise % {need to include 'sphere','cone','cylinder', but do not want to at this time since just tyring to get aspheres (10/21/2015)
                                    error('readLensFile:badLensType','Only spherical, aspherical, and plane surfaces can be handled with this funciton right now.')
                            end
                    end
                end
            end

    end
end

el.Orientation = eye(3);
el.Origin = zeros(3,1);
drawElementFirstTime(el);
end

function el = assignSurfProp(el,surfName,dat,datName, surfNum)
try
    el.Surfaces(surfNum).(surfName) = dat.surfaces{surfNum}.(lower(datName));
catch
    error('readLensFile:unknonwSphereParam','Missing parameter, %s, in surface %d.', datName, surfNum)
end
end

%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
