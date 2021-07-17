function Vecs = sphere_interp(V1, V2, samp_n)% So called SLERP
nV1 = V1 / norm(V1);
nV2 = V2 / norm(V2);
cosang = nV1 * nV2';
ang = acos(cosang);
if size(samp_n) == 1
    lingrid = linspace(0,1,samp_n)';
else
    lingrid = reshape(samp_n, [], 1);
end
Vecs = (sin((1-lingrid) * ang) * V1 + sin(lingrid * ang) * V2) / sin(ang);
end