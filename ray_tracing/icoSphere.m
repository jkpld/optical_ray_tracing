function [v,f] = icoSphere(n)
% create the vertices and faces to create a sphere out of trianbles. The
% number of base triangle sis 20, specifiing n = 2 will break of each of
% theses 20 triangles into four other triangles. 

% http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html

t = (1 + sqrt(5))/2;

v(12,3) = 0;
f(20,3) = 0;

% vertices
v(1,:) = [-1  t 0];
v(2,:) = [ 1  t 0];
v(3,:) = [-1 -t 0];
v(4,:) = [ 1 -t 0];

v(5,:) = [0 -1  t];
v(6,:) = [0  1  t];
v(7,:) = [0 -1 -t];
v(8,:) = [0  1 -t];

v(9,:)  = [ t 0 -1];
v(10,:) = [ t 0  1];
v(11,:) = [-t 0 -1];
v(12,:) = [-t 0  1];

v = bsxfun(@rdivide, v, sqrt(sum(v.*v,2))); % normalize the values.


% faces
f(1,:) = [0 11 5];
f(2,:) = [0 5 1];
f(3,:) = [0 1 7];
f(4,:) = [0 7 10];
f(5,:) = [0 10 11];

f(6,:) = [1 5 9];
f(7,:) = [5 11 4];
f(8,:) = [11 10 2];
f(9,:) = [10 7 6];
f(10,:) = [7 1 8];

f(11,:) = [3 9 4];
f(12,:) = [3 4 2];
f(13,:) = [3 2 6];
f(14,:) = [3 6 8];
f(15,:) = [3 8 9];

f(16,:) = [4 9 5];
f(17,:) = [2 4 11];
f(18,:) = [6 2 10];
f(19,:) = [8 6 7];
f(20,:) = [9 8 1];
f = f+1;

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
    ind = ind + 12;
end

end

%-%
%-% This is love: not that we loved God, but that he loved us and sent his
%-% Son as an atoning sacrifice for our sins. (1 John 4:10)
%-%
