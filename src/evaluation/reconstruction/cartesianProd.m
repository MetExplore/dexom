% From https://stackoverflow.com/questions/9834254/cartesian-product-in-matlab
function C = cartesianProd(varargin)
    args = varargin;
    n = nargin;

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end

    C = unique(G , 'rows');
end