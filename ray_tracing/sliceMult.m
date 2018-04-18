function out = sliceMult(a,b,c,gpu,transC)
% a b and c must all be the same size, the first and second dimensions must
% also be the same size, ex. size(a) = [n,n,m], size(a) = size(b) = size(c)
% out(:,:,i) = a(:,:,i)*b(:,:,i)*c(:,:,i)
%
% If c is empty, then just a(:,:,i)*b(:,:,i) will be calculated
%
% gpu is a boolian flag; if true, then the gpu will be used (if possible)
% to do the computation.

% if the 5th input is true, than c will be trasposed.

if nargin == 5
    trans = logical(transC);
else
    trans = false;
end

if nargin < 4
    error('sliceMult:badInput', 'must be four inputs to sliceMult')
end

if isempty(a) || isempty(b) 
    error('sliceMult:badInput', 'a and b cannot be empty')
end

if ~isempty(c)
    siz3 = [size(a,3),size(b,3),size(c,3)];
    if ~(numel(unique([size(a,1) size(a,2) size(b,1) size(b,2) size(c,1) size(c,2)])) == 1)
        error('sliceMult:badInput','The first 2 dimensions of a, b, and c must all have the same size.')
    end
else
    siz3 = [size(a,3),size(b,3)];
    if ~(numel(unique([size(a,1) size(a,2) size(b,1) size(b,2)])) == 1)
        error('sliceMult:badInput','The first 2 dimensions of a and b must all have the same size.')
    end
end

siz3u = unique(siz3);
if ~(numel(siz3u)==1)
    if numel(siz3u) > 2
        error('sliceMult:badInput', 'matrices a, b, and c must be the same size or some may have a 3rd dimension size of 1 (c can also be empty).')
    elseif ~any(siz3u==1)
        error('sliceMult:badInput', 'matrices a, b, and c must be the same size or some may have a 3rd dimension size of 1 (c can also be empty).')
    end
end


gpuError = 0;
out = [];
if gpu
    try
        if gpuDeviceCount > 0
            a_d = gpuArray(a);
            b_d = gpuArray(b);
                
            if isempty(c)
                out_d = pagefun(@mtimes,a_d,b_d);
            else
                c_d = gpuArray(c);
                if trans
                    c_d = pagefun(@ctranspose, c_d);
                    out_d = pagefun(@mtimes,b_d,c_d);
                else
                    out_d = pagefun(@mtimes,b_d,c_d);
                end
                
                out_d = pagefun(@mtimes,a_d,out_d);
            end

            wait(gpuDevice)
            out = gather(out_d);
            clear a_d b_d c_d out_d
        end
    catch
        gpuDevice([]);
        gpuError = 1;
        warning('sliceMult:gpuError','There was an error trying the calculation on the gpu. Calculation will be done on the cpu (if the computer didn''t just crash).');
    end
end

if ~gpuError && ~gpu
    
    
    
    out = zeros(size(b));

    if isempty(c)
        if all(siz3 == 1)
            out = a*b;
        elseif siz3(1) == 1 && siz3(2) ~= 1
            out = zeros(size(b));
            for i = 1:siz3(2)
                out(:,:,i) = a*b(:,:,i);
            end
        elseif siz3(1) ~= 1 && siz3(2) == 1
            out = zeros(size(a));
            for i = 1:siz3(1)
                out(:,:,i) = a(:,:,i)*b;
            end
        else
            out = zeros(size(a));
            for i = 1:siz3(1)
                out(:,:,i) = a(:,:,i)*b(:,:,i);
            end
        end
    else
        if trans
            if all(siz3 == 1)
                out = a*b*c';
            elseif siz3(1) == 1 && siz3(2) ~= 1 && siz3(3) ~= 1
                out = zeros(size(b));
                for i = 1:siz3(2)
                    out(:,:,i) = a*b(:,:,i)*c(:,:,i)';
                end
            elseif siz3(1) ~= 1 && siz3(2) == 1 && siz3(3) ~= 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b*c(:,:,i)';
                end
            elseif siz3(1) ~= 1 && siz3(2) ~= 1 && siz3(3) == 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b(:,:,i)*c';
                end
            elseif siz3(1) == 1 && siz3(2) == 1 && siz3(3) ~= 1
                out = zeros(size(c));
                for i = 1:siz3(3)
                    out(:,:,i) = a*b*c(:,:,i)';
                end
            elseif siz3(1) ~= 1 && siz3(2) == 1 && siz3(3) == 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b*c';
                end
            elseif siz3(1) == 1 && siz3(2) ~= 1 && siz3(3) == 1
                out = zeros(size(b));
                for i = 1:siz3(2)
                    out(:,:,i) = a*b(:,:,i)*c';
                end
            else
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b(:,:,i)*c(:,:,i)';
                end
            end
        else
            if all(siz3 == 1)
                out = a*b*c;
            elseif siz3(1) == 1 && siz3(2) ~= 1 && siz3(3) ~= 1
                out = zeros(size(b));
                for i = 1:siz3(2)
                    out(:,:,i) = a*b(:,:,i)*c(:,:,i);
                end
            elseif siz3(1) ~= 1 && siz3(2) == 1 && siz3(3) ~= 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b*c(:,:,i);
                end
            elseif siz3(1) ~= 1 && siz3(2) ~= 1 && siz3(3) == 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b(:,:,i)*c;
                end
            elseif siz3(1) == 1 && siz3(2) == 1 && siz3(3) ~= 1
                out = zeros(size(c));
                for i = 1:siz3(3)
                    out(:,:,i) = a*b*c(:,:,i);
                end
            elseif siz3(1) ~= 1 && siz3(2) == 1 && siz3(3) == 1
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b*c;
                end
            elseif siz3(1) == 1 && siz3(2) ~= 1 && siz3(3) == 1
                out = zeros(size(b));
                for i = 1:siz3(2)
                    out(:,:,i) = a*b(:,:,i)*c;
                end
            else
                out = zeros(size(a));
                for i = 1:siz3(1)
                    out(:,:,i) = a(:,:,i)*b(:,:,i)*c(:,:,i);
                end
            end
        end
    end
end


end
%-%
%-% For God so loved the world that he gave his one and only Son, that
%-% whoever believes in him shall not perish but have eternal life. (John
%-% 3:16)
%-%
