function [v,f] = icoCircle(n)
% Same as icoSphere but for a circle

% http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html


v(4,2) = 0;
f(1,3) = 0;


% v(1,:) = [0 0];
% v(2,:) = [1 0];
% v(3,:) = [0.5 sqrt(3)/2];
% 
% 
% v(:,1) = v(:,1) - 0.5;
% v(:,2) = v(:,2) - sqrt(3)/4;

% v(1,:) = [0 0];
v(1,:) = [0 1];
v(2,:) = ([cosd(120), -sind(120); sind(120), cosd(120)]*v(2,:)')';
v(3,:) = ([cosd(120), -sind(120); sind(120), cosd(120)]*v(3,:)')';

% v = bsxfun(@rdivide, v, sqrt(sum(v.*v,2))); % normalize the values.

f(1,:) = [1 2 3];
% f(2,:) = [1 3 4];
% f(3,:) = [1 4 2];

key = [];

for i = 1:(n-1)
   
    newfaces = zeros(size(f,1)*4,3);    
%     newfaces = [];
    for j = 1:size(f,1)
        [a, v, key] = getMiddlePoint(f(j,1), f(j,2),v,key);
        [b, v, key] = getMiddlePoint(f(j,2), f(j,3),v,key);
        [c, v, key] = getMiddlePoint(f(j,3), f(j,1),v,key);
        
        newfaces(1+4*(j-1),:) = [f(j,1), a, c];
        newfaces(2+4*(j-1),:) = [f(j,2), b, a];
        newfaces(3+4*(j-1),:) = [f(j,3), c, b];
        newfaces(4+4*(j-1),:) = [a, b, c];

    end
    
    f = newfaces;
    
end

end

function [ind, v, key] = getMiddlePoint(p1,p2,v,key)

smallerIdx = int64(min(p1,p2));
maxIdx = int64(max(p1,p2));

ckey = bitshift(smallerIdx, 32) + maxIdx;

ind = find(ismember(key,ckey));

if isempty(ind)
    pm = (v(p1,:)+v(p2,:))/2;
    pm = pm/sqrt(sum(pm.*pm));
    
    key(end+1) = ckey;
    
    v(end+1,:) = pm;
    
    ind = size(v,1);
else
    ind = ind + 3;
end

end

%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
