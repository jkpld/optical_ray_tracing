function [rayDat_out, surfNorms] = distanceToSphere(rayDat, R, sphereCent, n)
% distanceToSphere  propogate rays to sphere and calculate the sphere
% normals at the ray intersection points
%
% This function will calculate the distance along a ray that a ray needs to
% travel to intersect with a sphere of radius R centered at sphereCent. The
% distance to the first intersection point will be given. If there is no
% intersection point, then NaN will be used.
% This function will additionally output the surface normal vectors at each
% point of intersection as a n x 3 matrix in surfNorms.
%
% rayDat  -- n x 11 matrix of ray data, where n is the number of rays. The
%            11 columns are 
%            x,y,z,theta,phi,uv_x,uv_y,uv_z,w,l,I
%
% x,y,z -- starting position of ray
% theta -- polar angle of ray [-pi/2, pi/2]
% phi   -- azimuthal angle of ray [0, 2*pi)
% uv_x,uv_y,uv_z -- unit vector giving the direction of the ray
% w -- wavelength of ray
% l -- total optical path length the ray has traveled
% I -- amplitude of ray (the square of which gives the light intensity)
%
% R -- radius of the sphere
%
% sphereCent -- 3 x 1 vector giving the center of the sphere
%
% n -- refractive index of the material 

y = bsxfun(@minus, rayDat(:,1:3),sphereCent');

y_dot_uv = sum(y.*rayDat(:,6:8),2);

underSqrt = y_dot_uv.^2 - (sum(y.^2,2) - R^2);

good = underSqrt>=0;

h = NaN(size(rayDat,1),1);

h(good) = -y_dot_uv(good) - sign(R)*sqrt(underSqrt(good));

rayDat_out = rayDat;
rayDat_out(:,10) = h.*n;
rayDat_out(:,1:3) = rayDat(:,1:3) + bsxfun(@times, rayDat(:,6:8), h);

surfNorms = bsxfun(@minus, rayDat_out(:,1:3), sphereCent');
surfNorms = bsxfun(@rdivide, surfNorms, sqrt(sum(surfNorms.^2,2)));

end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
