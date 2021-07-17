function Vecs = SLERP(V1, V2, samp_t)% So called SLERP
nV1 = V1 / norm(V1);
nV2 = V2 / norm(V2);
cosang = nV1 * nV2';
ang = acos(cosang);
lingrid = reshape(samp_t, [], 1);
Vecs = (sin((1-lingrid) * ang) * V1 + sin(lingrid * ang) * V2) / sin(ang);
end