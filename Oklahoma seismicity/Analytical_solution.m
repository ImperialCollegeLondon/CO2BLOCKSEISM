
function [PD] = Analytical_solution(r,R,R_ext)
    
if r <= R
    if R <= R_ext  % open flow conditions or closed but pressure front has not yet reached boundaries
        PD = log(R/r);
    else
        PD = log(R_ext/r) + 2/2.25*(R/R_ext)^2 -3/4 ;
    end
else
    PD = 0;
end
end