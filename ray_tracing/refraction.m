function [rayDat,PtA,PtE,PrA] = refraction(rayDat, surfNorms, n1, n2)
% This function will calculate the rays new directional unit vectors after
% refracting through a surface with surface norms, surfNorms. The polar and
% azimuthal angles will also be updated. n1 is the refractive index of the
% medium the rays start in, and n2 is the refractive index of the medium
% the rays end in.
% Before calling this function, the rays should be propogated to the
% surface and the surface normals at the intersections should be calculated
% (surfNorms).
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
% n1 -- refractive index of first medium
% n2 -- refractive index of second medium
%
% PtA -- Polarization matrix for the electric field amplitudes. This is a
%        3 x 3 x n array where n is the number or rays. Using these, the
%        final electric field vector for ray i is given by 
%        E_final = PtA(:,:,i)*E.
%
% PtE -- Polarization matrix for the ray energy. This is a 3 x 3 x n array
%        where n is the number or rays. (This is very similar to PtA;
%        however, the Fresnel amplitude transmission coefficients have been
%        multiplied by the array scalling factor (n2.*cos2)/(n1.*cos1). The
%        energy of ray i, will now be given by |P(:,:,i)*E|^2, where E is
%        the electric field vector.)

% James Kapaldo, Jan. 2016

% 2016, 01, 23. Fixed an error in the k2 calculation (forgot to divide cos1
% by n2).

k1 = rayDat(:,4:6);
n1 = n1(:);
n2 = n2(:);

cos1 = sum(k1.*surfNorms,2);
cos2 = sign(cos1).*sqrt(1-(n1./n2).^2.*(1-cos1.^2));

% Find rays that will undergo total internal reflection
tir = imag(cos2) ~= 0; 

k2 = bsxfun(@times, k1, (n1./n2)) - bsxfun(@times, surfNorms, n1.*cos1./n2 - cos2);
k2 = bsxfun(@rdivide, k2, sqrt(sum(k2.^2,2)));

% ==== Note!!
% I will not include total internal reflection because the rays would no
% longer follow the sequential ray tracing that I use. If you want to
% include total internal reflection, then you can uncomment the below two
% lines marked, and also uncomment the line at the end of the next section.

% k2(tir,:) = k1(tir,:) - 2*bsxfun(@times, sum(k1(tir,:).*surfNorms(tir,:),2), surfNorms(tir,:));% Use for total internal reflection
% k2(tir,:) = nan(sum(tir),3); % Use to ignore total internal reflection


%% Calculate the polarization ray tracing matrix
% This uses the Fresnel reflection and refraction equations to calculate
% the transmitance and reflectace in the local coordinate system, and the
% transforms into absolute coordinates.
%
% See, "Three-dimensional polarization ray-tracing calculus I: definition
% and diattenuation", doi:10.1364/AO.50.002866

s1 = cross(k1,k2);

% For rays at normal incidance, we first need to generate a new s1 vector.
% I will just generate an arbitrary perpendicular vector using gram-schmidt
% orthogonalization.
s1length = sqrt(sum(s1.^2,2));
normInc = s1length == 0 & abs(cos1) == 1; % dont include tangent cases -- the cosine part

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
% cos2(tir) = cos1(tir);
ts = 2*n1.*cos1./(n1.*cos1 + n2.*cos2);
tp = 2*n1.*cos1./(n1.*cos2 + n2.*cos1);

rs = (n1.*cos1 - n2.*cos2)./(n1.*cos1+n2.*cos2);
rp = (n2.*cos1 - n1.*cos2)./(n2.*cos1+n1.*cos2);

% rs = @(n1,n2,th) (n1*cosd(th) - n2*sqrt(1-(n1/n2 * sind(th)).^2))./(n1*cosd(th)+n2*sqrt(1-(n1/n2 * sind(th)).^2));
% rp = @(n1,n2,th) (n2*cosd(th) - n1*sqrt(1-(n1/n2 * sind(th)).^2))./(n2*cosd(th)+n1*sqrt(1-(n1/n2 * sind(th)).^2));

JtA = zeros(3,3,n);
JtA(3,3,:) = 1;
JtA(1,1,:) = ts;
JtA(2,2,:) = tp;

JtE = zeros(3,3,n);
JtE(3,3,:) = 1;
JtE(1,1,:) = ts .* sqrt((n2.*cos2)./(n1.*cos1));
JtE(2,2,:) = tp .* sqrt((n2.*cos2)./(n1.*cos1));

JrA = zeros(3,3,n);
JrA(3,3,:) = 1;
JrA(1,1,:) = rs;
JrA(2,2,:) = rp;

PtA = sliceMult(Oout,JtA,inv_Oin,0);
PtE = sliceMult(Oout,JtE,inv_Oin,0);
PrA = sliceMult(Oout,JrA,inv_Oin,0);
% =========   Code for testing that P*k1 = k2
% test = zeros(n,6);
% testr = zeros(n,6);
% for i = 1:n
%     test(i,:) = [(PtA(:,:,i)*k1(i,:)')',k2(i,:)];
%     testr(i,:) = [(PtE(:,:,i)*k1(i,:)')',k2(i,:)];
% end

% Correction for tangent rays
PtA(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));
PtE(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));
PrA(:,:,cos1==0) = repmat(eye(3),1,1,sum(cos1==0));


% ==== Note!!
% I will not include total internal reflection because the rays would no
% longer follow the sequential ray tracing that I use. If you want to
% include total internal reflection, then you can switch the commented
% lines below; however, you must also add in the phase change of the
% differnent polarizations in the section of code below, which has not yet
% been programed.

k2(tir,:) = k1(tir,:) - 2 * bsxfun(@times, surfNorms(tir,:), sum(k1(tir,:).*surfNorms(tir,:),2)); % Use for total internal reflection
% k2(tir,:) = nan(sum(tir),3); % Use to ignore total internal reflection

rayDat(:,4:6) = k2;



end





%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
