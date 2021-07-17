function norms = norm_axis(Vecs, varargin)
if nargin < 2
    dim = 1 ;
else
    dim = varargin{1};
end
norms = sqrt(sum(Vecs.^2,dim));
end
