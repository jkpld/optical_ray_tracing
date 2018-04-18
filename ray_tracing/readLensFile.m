function el = readLensFile(file)
el = opticalElement3D;
el.Focal = zeros(1,2);
fid = fopen(file);
try
    if fid == -1
        error('readLensFile:badInput','Unable to open file, %s', file);
    else
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if strncmp(tline,'//',2)
                continue
            else
                
                [param,~,~,idx] = sscanf(tline,'%s:');

                switch param
                    case 'Manufacturer:'
                        el.Manufacture = tline(idx+1:end);
                    case 'ID:'
                        el.ElementID = tline(idx+1:end);
                    case 'Type:'
                        el.Type = lower(tline(idx+1:end));
                    case 'SubType:'
                        el.SubType = lower(tline(idx+1:end));
                    case 'crossSection:'
                        temp = lower(tline(idx+1:end));
                        if ~ismember({'circle','rectangle'},temp)
                            error('readLensFile:badCrossSection','The crossSection must be either "circle" or "rectangle".');
                        end
                        el.crossSection = temp; 
                    case 'crossSectionSize:'
                        temp =  sscanf(tline(idx+1:end),'%f');

                        if numel(temp)>2 || numel(temp)<1 || ~all(isfinite(temp))
                            error('readLensFile:badCrossSectionSize','The crossSectionSize must have 1 or 2 finite elements.');
                        end
            
                        if strcmp(el.crossSection,'circle')
                            temp = temp/2;
                        end
                        
                        el.crossSectionSize = temp(:)'; 
                    case 'sequential:'
                        temp =  sscanf(tline(idx+1:end),'%f');

                        if ~(numel(temp)==1 && (temp ==0 || temp ==1))
                            error('readLensFile:badSequential','Bad sequential, must be 0 or 1');
                        end

                        el.sequentialRayTrace = logical(temp);
                    case 'clearAperture:'
                        temp =  sscanf(tline(idx+1:end),'%f');

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
                    case 'f:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        if ~isempty(temp)
                            if numel(temp) ~= 1 || ~isfinite(temp)
                                error('readLensFile:badFocalLength','Bad focal length');
                            end
                        end
                        el.Focal(1) = temp;
                    case 'fb:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        if ~isempty(temp)
                            if numel(temp) ~= 1 || ~isfinite(temp)
                                error('readLensFile:badbackFocalLength','Bad back focal length');
                            end
                            el.Focal(2) = temp;
                        else
                            el.Focal(2) = [];
                        end
                                                
                    case 'NA:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        if ~isempty(temp)
                            if numel(temp) ~= 1 || ~isfinite(temp)
                                error('readLensFile:badNA','Bad NA value');
                            end
                        end
                        el.NA = temp;
                    case 'Price:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        el.Price = temp(1);
                    case 'DateOfInformation:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        el.DateOfInformation = temp;
                    case 'edgeThickness:'
                        temp =  sscanf(tline(idx+1:end),'%f');
                        if numel(temp) ~= 1 || ~isfinite(temp)
                            error('readLensFile:badedgeThickness','Bad edgeThickness');
                        end
                        edgeThickness = temp(1);
                    case 'numberOfSurfaces:'
                        numSurfaces =  sscanf(tline(idx+1:end),'%f');
                        if numel(numSurfaces) ~= 1
                            error('readLensFile:badFormat','The value of the parameter ''numberOfSurfaces'' must be a scalar number');
                        end
                        
                        el.Surfaces(numSurfaces) = nullSurface;
                        
                        for i = 1:numSurfaces
                            
                            % Header
                            tline = getNextLine(fid);
                            param = sscanf(tline,'%s');
                            if ~strcmp(param,'Surface')
                                error('readLensFile:badFormat','Bad file format in the surfaces sections, a surface must start with the line ''Surface''');
                            end
                            
                            % Material
                            tline = getNextLine(fid);
                            [param,~,~,idx] = sscanf(tline,'%s:');
                            if ~strcmp(param,'material:')
                                error('readLensFile:badFormat','Bad file format in the surfaces sections, Material must be the first parameter of a surface');
                            end
                            
                            if i ~= numSurfaces
                                el.Material{length(el.Material)+1} = sscanf(tline(idx+1:end),'%s');
                            end
                            
                            % surfaceType
                            tline = getNextLine(fid);
                            [param,~,~,idx] = sscanf(tline,'%s:');
                            if ~strcmp(param,'surfaceType:')
                                error('readLensFile:badFormat','Bad file format in the surfaces sections, zPlane must be the second parameter of a surface');
                            end
                            surfType = sscanf(tline(idx+1:end),'%s');
                            switch surfType
                                case 'sphere'
                                    el.Surfaces(i) = sphericalSurface;
                                    el.Surfaces(i).clearAperture = clearAperture;
                                    % R, radius
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    value = sscanf(tline(idx+1:end),'%f');
                                    switch param
                                        case 'R:'
                                            el.Surfaces(i).R = value;
                                        otherwise
                                            error('readLensFile:unknonwSphereParam','Unknown sphere parameter, %s',param)
                                    end
                                    
                                    % apex, 
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    value = sscanf(tline(idx+1:end),'%f');
                                    switch param
                                        case 'apex:'
                                            validateattributes(value,{'double'},{'finite','numel',3})
                                            el.Surfaces(i).Center = value(:) - [el.Surfaces(i).R;0;0];
                                        otherwise
                                            error('readLensFile:unknonwSphereParam','Unknown sphericalSurface parameter, %s',param)
                                    end
                                case 'plane'
                                    el.Surfaces(i) = planeSurface;
                                    el.Surfaces(i).clearAperture = clearAperture;
                                    for j = 1:2 
                                        tline = getNextLine(fid);
                                        [param,~,~,idx] = sscanf(tline,'%s:');
                                        value = sscanf(tline(idx+1:end),'%f');
                                        switch param
                                            case 'normVec:'
                                                validateattributes(value,{'double'},{'finite','numel',3})
                                                el.Surfaces(i).normVec = value(:);
                                            case 'point:'
                                                validateattributes(value,{'double'},{'finite','numel',3})
                                                el.Surfaces(i).Point = value(:);
                                            otherwise
                                                error('readLensFile:unknonwPlaneParam','Unknown planeSurface parameter, %s',param)
                                        end
                                    end
                                case 'backSurface'
                                    el.Surfaces(i) = backSurface;
                                    el.Surfaces(i).clearAperture = [];
                                    for j = 1:2 
                                        tline = getNextLine(fid);
                                        [param,~,~,idx] = sscanf(tline,'%s:');
                                        value = sscanf(tline(idx+1:end),'%f');
                                        switch param
                                            case 'normVec:'
                                                validateattributes(value,{'double'},{'finite','numel',3})
                                                el.Surfaces(i).normVec = value(:);
                                            case 'point:'
                                                validateattributes(value,{'double'},{'finite','numel',3})
                                                el.Surfaces(i).Point = value(:);
                                            otherwise
                                                error('readLensFile:unknonwBackSurfaceParam','Unknown backSurface parameter, %s',param)
                                        end
                                    end
                                case 'parabola'
                                    el.Surfaces(i) = parabolicSurface;
                                    el.Surfaces(i).clearAperture = clearAperture;
                                    good = false(1,4);
                                    for j = 1:4 
                                        tline = getNextLine(fid);
                                        [param,~,~,idx] = sscanf(tline,'%s:');
                                        value = sscanf(tline(idx+1:end),'%f');
                                        switch param
                                            case 'PFL:'
                                                el.Surfaces(i).PFL = value(:);
                                                good(1) = true;
                                            case 'EFL:'
                                                el.Surfaces(i).EFL = value(:);
                                                good(2) = true;
                                            case 'Angle:'
                                                el.Surfaces(i).Angle = value(:);
                                                good(3) = true;
                                            case 'baseDist:'
                                                el.Surfaces(i).baseDist = value(:);
                                                good(4) = true;
                                            otherwise
                                                error('readLensFile:unknonwParabolicParam','Unknown parabolicSurface parameter, %s',param)
                                        end
                                    end
                                    if ~all(good)
                                        error('readLensFile:badParabolaSurface','A reequired parameter for an Aperture surface is missing, must have "PFL","EFL","Angle","baseDist".');
                                    end
                                    
                                case 'aperture'
                                    
                                    el.Surfaces(i) = apertureSurface;
                                    el.Surfaces(i).clearAperture = clearAperture;
                                    el.Material(end) = [];
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    if ~strcmp(param,'shape:')
                                        error('readLensFile:badSurface','The first parameter of an aperture surface must be "shape".');
                                    else
                                        temp = lower(tline(idx+1:end));
                                        if ~ismember({'circle','rectangle'},temp)
                                            error('readLensFile:badApertureShape','The aperture shape must be either "circle" or "rectangle".');
                                        end
                                        el.Surfaces(i).Shape = temp;
                                    end
                        
                                    good = false(3,1);
                                    for j = 1:3 
                                        tline = getNextLine(fid);
                                        [param,~,~,idx] = sscanf(tline,'%s:');
                                        value = sscanf(tline(idx+1:end),'%f');
                                        switch param
                                            case 'normVec:'
                                                el.Surfaces(i).normVec = value(:);
                                                good(1) = true;
                                            case 'center:'
                                                el.Surfaces(i).Center = value(:);
                                                good(2) = true;
                                            case 'Extent:'
                                                el.Surfaces(i).Extent = value(:);
                                                good(3) = true;
                                            otherwise
                                                error('readLensFile:unknonwApertureParam','Unknown aperture parameter, %s',param)
                                        end
                                    end
                                    if ~all(good)
                                        error('readLensFile:badApertureSurface','A reequired parameter for an Aperture surface is missing, must have "shape","center","normVec","Extent".');
                                    end
                                case 'asphere'
                                    el.Surfaces(i) = asphericalSurface;
                                    
                                    % apex, 
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    value = sscanf(tline(idx+1:end),'%f');
                                    switch param
                                        case 'apex:'
                                            validateattributes(value,{'double'},{'finite','numel',3})
                                            el.Surfaces(i).Center = value(:);
                                        otherwise
                                            error('readLensFile:unknonwAsphereParam','Unknown asphericalSurface parameter, %s',param)
                                    end
                                    
                                    % R, radius
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    value = sscanf(tline(idx+1:end),'%f');
                                    switch param
                                        case 'R:'
                                            el.Surfaces(i).R = value;
                                        otherwise
                                            error('readLensFile:unknonwAsphereParam','Unknown asphericalSurface parameter, %s',param)
                                    end
                                    
                                    % k, radius
                                    tline = getNextLine(fid);
                                    [param,~,~,idx] = sscanf(tline,'%s:');
                                    value = sscanf(tline(idx+1:end),'%f');
                                    switch param
                                        case 'conic:'
                                            el.Surfaces(i).k = value;
                                        otherwise
                                            error('readLensFile:unknonwAsphereParam','Unknown asphericalSurface parameter, %s',param)
                                    end
                                    
                                    % A coefficients
                                    A = zeros(8,1);
                                    for j = 1:8
                                        tline = getNextLine(fid);
                                        [param,~,~,idx] = sscanf(tline,'%s:');
                                        value = sscanf(tline(idx+1:end),'%f');
                                        switch param
                                            case {'2nd:','4th:','6th:','8th:','10th:','12th:','14th:','16th:'}
                                                A(j) =value;
                                            otherwise
                                                error('readLensFile:unknonwAsphereParam','Unknown asphericalSurface parameter, %s',param)
                                        end
                                    end
                                    el.Surfaces(i).A = A;
                                    
                                otherwise % {need to include 'sphere','cone','cylinder', but do not want to at this time since just tyring to get aspheres (10/21/2015)
                                    error('readLensFile:badLensType','Only spherical, aspherical, and plane surfaces can be handled with this funciton right now.')
                            end
                            
                            % End Header
                            tline = getNextLine(fid);
                            if ~strcmp(tline,'end Surface')
                                error('readLensFile:badFormat','Bad file format in the surfaces sections, a surface must end with the line ''end Surface''');
                            end
                        end
                        
                    otherwise
                        error('readLensFile:unknownParameter','There was an unknown parameter, %s, in the lens file, %s.',param,file);
                end
            end
        end


    end
    

    
catch ME
    fclose(fid);
    throw(ME);
end

fclose(fid);

el.Orientation = eye(3);
el.Origin = zeros(3,1);
drawElementFirstTime(el);

end

function tline = getNextLine(fid)
    while 1
        tline = fgetl(fid);
        if ~ischar(tline)
            error('readLensFile:badFile','badFile format in the surfaces sections')
        end

        if strncmp(tline,'//',2)
            continue
        else
            break;
        end
    end
end

%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
