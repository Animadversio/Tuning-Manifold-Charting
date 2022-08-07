function y = sem(varargin)
% function y = sem(x,dim)
% Standard error of the mean, where x is a matrix or vector
% dim: dimension of operation
% Should handle nans well.
if nargin == 1
    x = varargin{1};
    [~,dim] = max(size(x)); 
elseif nargin == 2
    x = varargin{1};
    dim = varargin{2};
end

nSamples_notNan = sum(~isnan(x),dim)-1;

y = nanstd(x,[],dim) ./ sqrt(nSamples_notNan);

end