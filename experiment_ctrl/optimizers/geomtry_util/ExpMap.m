function y = ExpMap(x, tang_vec)
	% Modified Feb. 25th to allow x of arbitrary norm.
     EPS = 1E-10;
%     assert(abs(norm(x)-1) < EPS);
    assert(abs(norm(x))>EPS)
    unit_x = x ./ norm(x);
    % assert(sum(x * tang_vec') < EPS);
    angle_dist = sqrt(sum(tang_vec.^2, 2)); % vectorized
    uni_tang_vec = tang_vec ./ angle_dist;
    unit_x = repmat(unit_x, size(tang_vec, 1), 1); % vectorized
    y = cos(angle_dist) .* unit_x + sin(angle_dist) .* uni_tang_vec;
    y = y * norm(x);
end