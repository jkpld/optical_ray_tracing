function [rayDat,P] = reflection(rayDat, surfNorms,n1,n2)
% This function will calculate the rays new direction unit vectors after
% reflecting off a surface with surface norms, surfNorms. The polar and
% azimuthal angles will also be updated.
% Before calling this function, the rays should be propogated to the
% surface and the surface normals at the intersections should be calculated
% (surfNorms).
% If n2 or both n1 and n2 are empty arrays ([]), then the material being
% reflected off is assumed to be a perfect conductor
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
% surfNorms -- n x 3 matrix of the surface normals at the points of
% intersection. n is the number of rays = size(rayDat,1)
%
% n1 -- refractive index of the medium the ray is traveling in (ex. air,
%       glass)
%
% n2 -- refractive index of the medium the rays are reflecting off (ex.
%       metal)
%
% P -- Polarization matrix for the electric field amplitudes. This is a
%      3 x 3 x n array where n is the number or rays. Using these, the
%      final electric field vector for ray i is given by 
%      E_final = P(:,:,i)*E. Also, since this is reflection, the intensity
%      of the rays are just given by |E_final|^2.
%
% Note! Reflection is assumed to be from perfect conductors 

k1 = rayDat(:,4:6);
k2 = k1 - 2 * repmat(sum(k1.*surfNorms,2),1,3) .* surfNorms;

rayDat(:,4:6) = k2;


%% Calculate the polarization ray tracing matrix
% This uses the Fresnel reflection and refraction equations to calculate
% the transmitance and reflectace in the local coordinate system, and the
% transforms into absolute coordinates.
%
% See, "Three-dimensional polarization ray-tracing calculus I: definition
% and diattenuation", doi:10.1364/AO.50.002866

cos1 = sum(k1.*surfNorms,2);
cos2 = sign(cos1).*sqrt(1-(n1/n2)^2*(1-cos1.^2));

k2 = (n1/n2)*k1 + surfNorms.*repmat((cos2 - n1.*cos1),1,3);
k2 = k2./repmat(sqrt(sum(k2.^2,2)),1,3);

s1 = cross(k1,k2);

% For rays at normal incidance, we first need to generate a new s1 vector.
% I will just generate an arbitrary perpendicular vector using gram-schmidt
% orthogonalization.
length = sqrt(sum(s1.^2,2));
normInc = length == 0 & abs(cos1) == 1; % dont include tangent cases -- the cosine part

s1(normInc,:) = gramSchmidt1(rand(sum(normInc),3),k1(normInc,:));

s1 = s1./repmat(sqrt(sum(s1.^2,2)),1,3);

p1 = cross(k1,s1);
p2 = cross(k2,s1);

% s2 = s1;

% Store the coordinate transformation matrices as 3x3xn matrix where n is
% the number of rays

n = size(s1,1);

inv_Oin(3,3,:) = k1(:,3);
inv_Oin(3,2,:) = k1(:,2);
inv_Oin(3,1,:) = k1(:,1);
inv_Oin(2,3,:) = p1(:,3);
inv_Oin(2,2,:) = p1(:,2);
inv_Oin(2,1,:) = p1(:,1);
inv_Oin(1,3,:) = s1(:,3);
inv_Oin(1,2,:) = s1(:,2);
inv_Oin(1,1,:) = s1(:,1);

Oout(3,3,:) = k2(:,3);
Oout(2,3,:) = k2(:,2);
Oout(1,3,:) = k2(:,1);
Oout(3,2,:) = p2(:,3);
Oout(2,2,:) = p2(:,2);
Oout(1,2,:) = p2(:,1);
Oout(3,1,:) = s1(:,3);
Oout(2,1,:) = s1(:,2);
Oout(1,1,:) = s1(:,1);

% Calculate transmission coefficients

% ts = 2*n1*cos1./(n1*cos1 + n2*cos2);
% tp = 2*n1*cos1./(n1*cos2 + n2*cos1);

rs = (n1.*cos1 - n2.*cos2)./(n1.*cos1+n2.*cos2);
rp = (n2.*cos1 - n1.*cos2)./(n2.*cos1+n1.*cos2);

JrA = zeros(3,3,n);
JrA(3,3,:) = 1;
JrA(1,1,:) = rs;
JrA(2,2,:) = rp;

P = sliceMult(Oout,JrA,inv_Oin,0);

% =========   Code for testing that P*k1 = k2
% test = zeros(n,6);
% testr = zeros(n,6);
% for i = 1:n
%     test(i,:) = [(PtA(:,:,i)*k1(i,:)')',k2(i,:)];
%     testr(i,:) = [(PtE(:,:,i)*k1(i,:)')',k2(i,:)];
% end

% Correction for tangent rays
P(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));

end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
