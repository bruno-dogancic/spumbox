function x = dynlbl(a,n,type)

% Creates dynamical label (variable names) x with label name 'a'
% 'a' should be string of desired variable name
%           x = dynlbl(a,n,type)
% If n is a scalar then the resulting variable names are:
%
%     Type 1: a1, a2, a3 ... an
%     Type 2: a(1), a(2), a(3) ... a(n)
%
% If n is a vector [q,p] then the resulting names are:
%
%     Type 1: [a11, a12, a13 ... a1p;
%              a21, a22, a23 ... a2p;
%               |   ...  ... ... ...;
%              aq1, aq2, aq3 ... aqp],
%
%     Type 2: [a(1,1), a(1,2), a(1,3) ... a(1,p);
%              a(2,1), a(2,2), a(2,3) ... a(2,p);
%                |       ...    ...   ...   ... ;
%              a(q,1), a(q,2), a(q,3) ... a(q,p)],
%
%     Type 3: Works only if 'n' is a vector, and transforms 
%             Type 1 into row vector:
%             [a11, a12, a13 ... a1p, a21, a22, a23 ... a2p, 
%              aq1, aq2, aq3 ... aqp],
%
%     Type 4: Works only if 'n' is a vector, and transforms 
%             Type 2 into row vector:
%             [a(1,1), a(1,2), ... a(1,p), a(2,1), a(2,2), ... a(2,p),
%              a(q,1), a(q,2), ... a(q,p)].


if length(n) == 1
    x = cell(n,length(n));
    if type == 1
        
       for i = 1:n
            x{i} = sprintf('%s%d',a,i);
       end 
    elseif type == 2
        for i = 1:n
            x{i} = sprintf('%s(%d)',a,i);
        end
    else
        error('For scalar n, type should be either 1 or 2')
    end
elseif length(n) == 2
    x = cell(n(1),n(2));
    if type == 1
        for i = 1:n(1)
            for j = 1:n(2)
                x{i,j} = sprintf('%s%d%d',a,i,j);
            end
        end
    elseif type == 2
        for i = 1:n(1)
            for j = 1:n(2)
                x{i,j} = sprintf('%s(%d,%d)',a,i,j);
            end
        end
    elseif type == 3
        for i = 1:n(1)
            for j = 1:n(2)
                x{i,j} = sprintf('%s%d%d',a,i,j);
            end
        end
        x = transpose(reshape(transpose(x),1,[])); % convert matrix to column vector by first horz cat each row
    elseif type == 4
        for i = 1:n(1)
            for j = 1:n(2)
                x{i,j} = sprintf('%s(%d,%d)',a,i,j);
            end
        end
        x = transpose(reshape(transpose(x),1,[])); % convert matrix to column vector by first horz cat each row   
    else
        error('For vector n, type should be either 1, 2, 3 or 4')
    end
else
    error('n should be either scalar or vector of length 2')
end
