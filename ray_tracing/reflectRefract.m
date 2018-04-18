function [rayDat,PA,PE] = reflectRefract(rayDat, surfNorms, n1, n2)
% reflectRefract   Reflect or refract a ray off a surface.
%
% This function will either reflect or refract a ray off a surface with
% surface normals, surfNorms. The ray undergoes reflection if the the
% unpolarized reflectance is greater than half, R > 0.5. Before calling this
% function, the rays should be propogated to the surface and the surface
% normals at the intersections should be calculated (surfNorms).
%
% rayDat  -- n x 8 matrix of ray data, where n is the number of rays. The
%            10 columns are 
%            x,y,z,uv_x,uv_y,uv_z,w,l
%
%               x,y,z -- starting position of ray
%               uv_x,uv_y,uv_z -- unit vector giving the direction of the
%                                 ray
%               w -- wavelength of ray
%               l -- total optical path length the ray has traveled 
% 
% surfNorms -- n x 3 matrix of the surface normals at the points of
%              intersection. n is the number of rays = size(rayDat,1)
%
% n1 -- refractive index of first medium
% n2 -- refractive index of second medium
%
% PA -- Polarization matrix for the electric field amplitudes. This is a
%       3 x 3 x n array where n is the number or rays. Using these, the
%       final electric field vector for ray i is given by 
%       E_final = PtA(:,:,i)*E.
%
% PE -- Polarization matrix for the ray energy. This is a 3 x 3 x n array
%       where n is the number or rays. (This is very similar to PtA;
%       however, the Fresnel amplitude transmission coefficients have been
%       multiplied by the array scalling factor (n2.*cos2)/(n1.*cos1). The
%       energy of ray i, will now be given by |P(:,:,i)*E|^2, where E is
%       the electric field vector.)

% James Kapaldo, 2016, 1, 24

k1 = rayDat(:,4:6);
n1 = n1(:);
n2 = n2(:);

% Refract all of the rays

cos1 = sum(k1.*surfNorms,2);
cos2 = sign(cos1).*sqrt(1-(n1./n2).^2.*(1-cos1.^2));

k2 = bsxfun(@times, k1, (n1./n2)) - bsxfun(@times, surfNorms, n1.*cos1./n2 - cos2);
k2 = bsxfun(@rdivide, k2, sqrt(sum(k2.^2,2)));


%% Calculate the polarization ray tracing matrix
% This uses the Fresnel reflection and refraction equations to calculate
% the transmitance and reflectace in the local coordinate system, and the
% transforms into absolute coordinates.
%
% See, "Three-dimensional polarization ray-tracing calculus I: definition
% and diattenuation", doi:10.1364/AO.50.002866

s1 = crossProd(k1,k2);

% For rays at normal incidance, we first need to generate a new s1 vector.
% I will just generate an arbitrary perpendicular vector using gram-schmidt
% orthogonalization.
s1_L = sqrt(sum(s1.*s1,2));
normInc = s1_L == 0 & abs(cos1) == 1; % dont include tangent cases -- the cosine part

s1(normInc,:) = gramSchmidt1(rand(sum(normInc),3),k1(normInc,:));

s1_L(normInc) = sqrt(sum(s1(normInc,:).*s1(normInc,:),2));

s1 = bsxfun(@rdivide, s1, s1_L);

p1 = crossProd(k1,s1);
p2 = crossProd(k2,s1);

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

ts = 2*n1.*cos1./(n1.*cos1 + n2.*cos2);
tp = 2*n1.*cos1./(n1.*cos2 + n2.*cos1);

rs = (n1.*cos1 - n2.*cos2)./(n1.*cos1 + n2.*cos2);
rp = (-n2.*cos1 + n1.*cos2)./(n2.*cos1 + n1.*cos2);

% rs = @(n1,n2,th) (n1*cosd(th) - n2*sqrt(1-(n1/n2 * sind(th)).^2))./(n1*cosd(th)+n2*sqrt(1-(n1/n2 * sind(th)).^2));
% rp = @(n1,n2,th) (n2*cosd(th) - n1*sqrt(1-(n1/n2 * sind(th)).^2))./(n2*cosd(th)+n1*sqrt(1-(n1/n2 * sind(th)).^2));

reflected = rs.*conj(rs) + rp.*conj(rp) > 1;

JA = zeros(3,3,n);
JA(3,3,:) = 1;

JA(1,1,~reflected) = ts(~reflected);
JA(2,2,~reflected) = tp(~reflected);
JA(1,1,reflected) = rs(reflected);
JA(2,2,reflected) = rp(reflected);

JE = JA;

JE(1,1,~reflected) = ts(~reflected) .* sqrt((n2(~reflected).*cos2(~reflected))./(n1(~reflected).*cos1(~reflected)));
JE(2,2,~reflected) = tp(~reflected) .* sqrt((n2(~reflected).*cos2(~reflected))./(n1(~reflected).*cos1(~reflected)));

% Coordinate transform the polarization matrices from the surface normal
% coordinate system to the element coordinate system.
PA = sliceMult(Oout,JA,inv_Oin,0);
PE = sliceMult(Oout,JE,inv_Oin,0);

% =========   Code for testing that P*k1 = k2
% test = zeros(n,6);
% testr = zeros(n,6);
% for i = 1:n
%     test(i,:) = [(PA(:,:,i)*k1(i,:)')',k2(i,:)];
%     testr(i,:) = [(PE(:,:,i)*k1(i,:)')',k2(i,:)];
% end
% test
% testr

% Correction for tangent rays
PA(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));
PE(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));

% Calculate reflected rays directions.

k2(reflected,:) = k1(reflected,:) - 2 * bsxfun(@times, surfNorms(reflected,:), sum(k1(reflected,:).*surfNorms(reflected,:),2)); 


rayDat(:,4:6) = k2;

end

function c = crossProd(a,b)
c = [a(:,2).*b(:,3) - a(:,3).*b(:,2), a(:,3).*b(:,1) - a(:,1).*b(:,3), a(:,1).*b(:,2) - a(:,2).*b(:,1)];
end
%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
