function Vecs = TangLERP(refV, tangV, samp_t)% So called SLERP
% EPS = 1E-6;
% norm_ref = norm(refV);
norm_tan = norm(tangV);
lingrid = reshape(samp_t, [], 1);
Vecs = refV  + lingrid * tangV / norm_tan;
end