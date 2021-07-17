function Vecs = TangSLERP(refV, tangV, samp_th)% So called SLERP
EPS = 1E-6;
norm_ref = norm(refV);
assert(abs(norm_ref)>EPS)
tangV_or = tangV - tangV * (tangV*refV') / norm_ref.^2;
norm_tan = norm(tangV_or);
lingrid = reshape(samp_th, [], 1);
Vecs = (sin((1-lingrid)) * refV / norm_ref + sin(lingrid) * tangV_or / norm_tan) * norm_ref;
end