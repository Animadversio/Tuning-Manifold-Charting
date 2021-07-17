function Vecs = LERP(V1, V2, samp_t)% So called SLERP
lingrid = reshape(samp_t, [], 1);
Vecs = ((1-lingrid) * V1 + lingrid * V2);
end