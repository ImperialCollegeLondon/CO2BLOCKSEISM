
function [PD] = Analytical_solution(r,R,R_ext)

if r <= R
    % open flow conditions or closed but pressure front has not yet reached boundaries
    if R <= R_ext
        PD = log(R/r);
    else
        PD = log(R_ext/r) + 2/2.25*(R/R_ext)^2 -3/4 ;
    end
else
    PD = 0;
end
end
