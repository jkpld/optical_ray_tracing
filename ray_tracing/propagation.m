function rayDat = propagation(rayDat, d, n)
% This function will propagate the rays a distance d through a medium with
% refractive index n. The refractive index is used to keep track of the
% optical path length traveled.
%
% rayDat  -- n x 8 matrix of ray data, where n is the number of rays. The
%            11 columns are 
%            x,y,z,uv_x,uv_y,uv_z,w,l
%
% x,y,z -- starting position of ray
% uv_x,uv_y,uv_z -- unit vector giving the direction of the ray
% w -- wavelength of ray
% l -- total optical path length the ray has traveled
%
% n -- refractive index of medium

if numel(d) == 1
    rayDat(:,1:3) = rayDat(:,1:3) + rayDat(:,4:6)*d;
else
    rayDat(:,1:3) = rayDat(:,1:3) + rayDat(:,4:6).*repmat(d,1,3);
end
rayDat(:,8) = rayDat(:,8) + n.*d;


end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
