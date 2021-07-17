function ang = ang_dist(V1,V2)
nV1 = V1 ./ norm_axis(V1,2);
nV2 = V2 ./ norm_axis(V2,2);
cosang = nV1 * nV2';
ang = real(acos(cosang));
end
