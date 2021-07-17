function [vtang_old, vtang_new] = InvExpMap(xold, xnew)
    xold = xold ./ norm(xold);
    xnew = xnew ./ norm(xnew);
    angle = acos(dot(xold, xnew)); % 0-pi angle distance between old and new position. 
    vtang_old = xnew - cos(angle) * xold; 
    vtang_old = vtang_old / norm(vtang_old) * angle;
    vtang_new = - xold + cos(angle) * xnew; 
    vtang_new = vtang_new / norm(vtang_new) * angle;
end