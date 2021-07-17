function codes_renorm = renormalize(codes, norms)
% norms should be a scaler or a 1d array same length as codes' first
% dimension
    assert(size(codes,1) == length(norms) || length(norms) == 1)
    norms = reshape(norms, [],1);
    code_norm = sqrt(sum(codes.^2, 2));
    codes_renorm = norms .* codes ./ code_norm;  % norms should be a 1d array
return 