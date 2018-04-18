function [intersects, inside, onEdge, dmin, dmax] = rayIntersectAABB(ray,box)

% http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection

% 2016, 01, 23


%First row of box in minimum point and second row is maximum point of
%an axis aligned box.

rayInvDir = 1./ray(:,4:6);
rayInvDirSign = (sign(rayInvDir) < 0) + 1;
% if ray direction is negative then then sign is 2, if ray direction is
% positive then the sign is 1

dmin = (box(rayInvDirSign(:,1),1) - ray(:,1)) .* rayInvDir(:,1);
dmax = (box(3-rayInvDirSign(:,1),1) - ray(:,1)) .* rayInvDir(:,1);

dymin = (box(rayInvDirSign(:,2),2) - ray(:,2)) .* rayInvDir(:,2);
dymax = (box(3-rayInvDirSign(:,2),2) - ray(:,2)) .* rayInvDir(:,2);

doesntIntersect = dmin > dymax | dmax < dymin;

dmin(~doesntIntersect) = max(dmin(~doesntIntersect),dymin(~doesntIntersect));
dmax(~doesntIntersect) = min(dmax(~doesntIntersect),dymax(~doesntIntersect));

dzmin = inf(size(doesntIntersect));
dzmax = -inf(size(doesntIntersect));

dzmin(~doesntIntersect,1) = (box(rayInvDirSign(~doesntIntersect,3),3) - ray(~doesntIntersect,3)) .* rayInvDir(~doesntIntersect,3);
dzmax(~doesntIntersect,1) = (box(3-rayInvDirSign(~doesntIntersect,3),3) - ray(~doesntIntersect,3)) .* rayInvDir(~doesntIntersect,3);

doesntIntersect = doesntIntersect | (dzmin > dmax) | (dzmax < dmin);

dmin(~doesntIntersect) = max(dmin(~doesntIntersect),dzmin(~doesntIntersect));
dmax(~doesntIntersect) = min(dmax(~doesntIntersect),dzmax(~doesntIntersect));

intersects = (~doesntIntersect) & (dmin > 0) & (dmax > 0);

onEdge = dmin == 0 | dmax == 0;

inside = xor(dmin < 0,dmax < 0) & ~onEdge;

end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
