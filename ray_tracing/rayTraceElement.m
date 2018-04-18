function [rays, PA, PE, terminated] = rayTraceElement(el,rayDat)


PA = [];
PE = [];
terminated = [];

if size(rayDat,3) > 1
    error('rayTraceElement:badInput','The input ray data should only have one page, size(rayDat,3) = 1.')
end

numRys = size(rayDat,1);

if numRys > 3000
    useGPU = 1;
else
    useGPU = 0;
end

numSurfs = 2*sum(~(isa(el.Surfaces,'apertureSurface') | isa(el.Surfaces,'backSurface'))) + sum(isa(el.Surfaces,'apertureSurface'));

rays = nan(numRys,8,numSurfs+1);
rays(:,:,1) = rayDat;
% Transform rays into coordinate system of lens
rays(:,:,1) = transformRays(el,rays(:,:,1));

% Determine if the rays intersect with the elements bounding box, which is
% now axis aligned since the rays have been put in the elements coordinate
% system.
[intersects, ~, ~, dmin, ~] = rayIntersectAABB(rays(:,:,1), [el.ExtentPrimary.Vertices(:,1)';el.ExtentPrimary.Vertices(:,8)']);
if ~any(intersects)
    terminated = true(size(rayDat,1),1);
    return;    
end


% tempRays = rays(intersects,1:3,1) + bsxfun(@times, rays(intersects,4:6,1), dmin(intersects));
% [~,inside] = removeRaysOutsideCrossSection(el,tempRays);
% 
% 
% intersects(intersects) = intersects(intersects) & inside;

% intersects = intersects | onEdge;
% Get the distance to the first and the last surface and use these distance
% to descide which side of the element the rays should enter from.

h1 = goToSurf(el.Surfaces(1), rays(intersects,:,1), 1);
h2 = goToSurf(el.Surfaces(end), rays(intersects,:,1), 1);

h1(isnan(h1)) = inf;
h2(isnan(h2)) = inf;

forward = false(numRys,1);
reverse = forward;

forward(intersects) = h1 < h2;
reverse(intersects) = ~forward(intersects);%h2 < h1;

terminated = ~intersects;

% Pre-calculate the index of refractions.
if any(forward) || any(reverse)
    n(size(rays,1),length(el.Material)+2) = 0;
    n(:,1) = 1;
    n(:,end) = 1;
    for i = 1:length(el.Material)
        n(:,i+1) = refractiveIndex(el.Material{i},rays(:,7,1))';
    end
end


% Allocate space for the polarization matrices
% numSurf = sum(~isa(el.Surfaces,'apertureSurface') | ~isa(el.Surfaces,'backSurface'));
PA = zeros(3,3,numRys);
PE = PA;

surfs = el.Surfaces;


% Ray trace all rays entering the first surface of the element
if any(forward)

    norms = zeros(length(forward),3);
    inside = false(length(forward),1);
    crntRys = forward;
    
    matIdx = 1;
    counter = 1;
    
    for i = 1:length(surfs)

        if any(crntRys)

            switch class(surfs(i))
                case 'apertureSurface'
                    
                    [~,rays(crntRys,:,counter+1),~, inside(crntRys)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el); 

%                     crntRys = forward & ~terminated;

                    counter = counter + 1;
                    
                case 'backSurface'
                    if i == 1
                        [~, rays(crntRys,:,counter+1)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el);
                        terminated = terminated & forward;
                        break;
                    else
                        continue;
                    end
                otherwise

                    [~,rays(crntRys,:,counter+1),N, inside(crntRys)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el); 

                    if ~isempty(N)
                        norms(crntRys,:) = N;
                    end
                    
                    terminated = terminated & ~inside;
                    crntRys = forward & ~terminated;

                    if any(crntRys)

                        [rays(crntRys,:,counter+2),PA_temp,PE_temp] = reflectRefract(rays(crntRys,:,counter+1),norms(crntRys,:),n(crntRys,matIdx),n(crntRys,matIdx+1));
                    else
                        break;
                    end
                    
                    if matIdx == 1
                        PA(:,:,crntRys) = PA_temp;
                        PE(:,:,crntRys) = PE_temp;
                    else
                        PA(:,:,crntRys) = sliceMult(PA_temp, PA(:,:,crntRys), [], useGPU);
                        PE(:,:,crntRys) = sliceMult(PE_temp, PE(:,:,crntRys), [], useGPU);
                    end
                    
                    matIdx = matIdx + 1;
                    counter = counter + 2;

            end

            terminated(crntRys) = terminated(crntRys) | any(isnan(rays(crntRys,:,counter)),2);
            crntRys = forward & ~terminated;
            
        else
            break;
        end

    end
end



if any(reverse)

    norms = zeros(length(reverse),3);
    inside = false(length(reverse),1);
    crntRys = reverse & ~terminated;

    surfs = flip(surfs);
    n = flip(n,2);
    
    matIdx = 1;
    counter = 1;
    
    for i = 1:length(surfs)

        if any(crntRys)
            
            
            switch class(surfs(i))
                case 'apertureSurface'
                    
                    [~,rays(crntRys,:,counter+1),~, inside(crntRys)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el); 
%                     crntRys = reverse & ~terminated;

                    counter = counter + 1;
%                     continue;
                case 'backSurface'
                    if i == 1
                        [~, rays(crntRys,:,counter+1)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el);
                        terminated = terminated & reverse;
                        break;
                    else
                        continue;
                    end
                otherwise

                    [~,rays(crntRys,:,counter+1),N, inside(crntRys)] = goToSurf(surfs(i), rays(crntRys,:,counter), n(crntRys,matIdx), el); 

                    if ~isempty(N)
                        norms(crntRys,:) = N;
                    end
                    
                    terminated = terminated & ~inside;
                    crntRys = reverse & ~terminated;

                    
                    
                    if any(crntRys)
                        [rays(crntRys,:,counter+2),PA_temp,PE_temp] = reflectRefract(rays(crntRys,:,counter+1),norms(crntRys,:),n(crntRys,matIdx),n(crntRys,matIdx+1));
                    else
                        break;
                    end
                    
                    if matIdx == 1
                        PA(:,:,crntRys) = PA_temp;
                        PE(:,:,crntRys) = PE_temp;
                    else
                        PA(:,:,crntRys) = sliceMult(PA_temp, PA(:,:,crntRys), [], useGPU);
                        PE(:,:,crntRys) = sliceMult(PE_temp, PE(:,:,crntRys), [], useGPU);
                    end
                    
                    matIdx = matIdx + 1;
                    counter = counter + 2;
                    

            end
                        
            terminated(crntRys) = terminated(crntRys) | any(isnan(rays(crntRys,:,counter)),2);
            crntRys = reverse & ~terminated;
                        
        else
            break;
        end
    end
end





rays(:,:,1) = []; % remove the first page, no need to output it since we will have it in the calling function as what we originally sent.

rays(:,:,:) = invTransformRays(el,rays(:,:,:));

if any(~terminated)
    PA(:,:,~terminated) = sliceMult(el.Orientation, PA(:,:,~terminated), el.Orientation', useGPU);
    PE(:,:,~terminated) = sliceMult(el.Orientation, PE(:,:,~terminated), el.Orientation', useGPU);
else
    PA = [];
    PE = [];
end
end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
