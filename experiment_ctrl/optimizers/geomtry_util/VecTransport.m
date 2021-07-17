function v_transport = VecTransport(xold, xnew, v)
    xold = xold ./ norm(xold);
    xnew = xnew ./ norm(xnew);
    x_symm_axis = xold + xnew;
    v_transport = v - 2 * dot(v, x_symm_axis)/norm(x_symm_axis)^2 * x_symm_axis; % Equation for vector parallel transport along geodesic
end