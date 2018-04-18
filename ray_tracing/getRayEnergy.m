function E = getRayEnergy(PE, rayDir)


numRys = size(rayDir,1);
if numRys > 3000
    useGPU = 1;
else
    useGPU = 0;
end

uv2 = gramSchmidt1(rand(size(rayDir,1),3),rayDir);
uv2 = bsxfun(@rdivide, uv2, sqrt(sum(uv2.*uv2,2)));

uv3 = crossProd(rayDir,uv2);

T(3,3,:) = uv3(:,3);
T(2,3,:) = uv3(:,2);
T(1,3,:) = uv3(:,1);
T(3,2,:) = uv2(:,3);
T(2,2,:) = uv2(:,2);
T(1,2,:) = uv2(:,1);
T(3,1,:) = rayDir(:,3,1);
T(2,1,:) = rayDir(:,2,1);
T(1,1,:) = rayDir(:,1,1);

% Tt = permute(T,[2,1,3]);
J = sliceMult(T,0.5*diag([0 1 1]), T,useGPU,true);

% PE_dag = conj(permute(PE,[2,1,3]));

Jf = sliceMult(PE,J,PE,useGPU,true);
E = squeeze(Jf(1,1,:) + Jf(2,2,:) + Jf(3,3,:));

end

function c = crossProd(a,b)
c = [a(:,2).*b(:,3) - a(:,3).*b(:,2), a(:,3).*b(:,1) - a(:,1).*b(:,3), a(:,1).*b(:,2) - a(:,2).*b(:,1)];
end
%-%
%-% But he was pierced for our transgressions, he was crushed for our
%-% iniquities; the punishment that brought us peace was on him, and by
%-% his wounds we are healed. (Isaiah 53:5)
%-%
