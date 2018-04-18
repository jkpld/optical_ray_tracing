function out = gramSchmidt1(v,u)
% Use Gram-Schmidt orthogonalization to make one vector perpendicular to
% another. Take in two vectors, and make the first perpendicular to the
% second.
%
% Note, the output vector is not normalized; it is not a unit vector.
%
% James Kapaldo, Jan. 2016

out = v - repmat(sum(v.*u,2)./sum(u.*u,2),1,3) .* u;

end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
