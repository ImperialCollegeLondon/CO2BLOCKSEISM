%#ok<*LLMNC>
function [Sigma_n,Tau] = stress_projection(delta,SH_max,Sh_min,Sv,Azi,Dip)

    L_PG = [cos(deg2rad(delta)) sin(deg2rad(delta)) 0;...         % delta is the rotation angle from SHmax to the N direction
            sin(deg2rad(delta)) -cos(deg2rad(delta)) 0;...        % Rotation matrix: principal stress (SH-Sh-Sv) to geographical system (N-E-D)
            0 0 -1] ;                                             % Assuming that Sv is one of the principal stresses

    Sigma = [SH_max 0 0; 0 Sh_min 0; 0 0 Sv] ;                    % principal stress tensor

    Sigma_G = L_PG * Sigma * L_PG';                               % stress tensor in the geographical coordinate

    n_n = [-(sin(deg2rad(Azi)))*(sin(deg2rad(Dip)));...           % normal of the fault (in N-E-d system)
           (cos(deg2rad(Azi)))*(sin(deg2rad(Dip)));...            % normal to the footwall
           -cos(deg2rad(Dip))];

    n_s = [cos(deg2rad(Azi)) ; sin(deg2rad(Azi)) ; 0];            % strike direction (in N-E-d system)

    n_d = [-(sin(deg2rad(Azi)))*(cos(deg2rad(Dip)));...           % dip direction (in N-E-d system)
           (cos(deg2rad(Azi)))*(cos(deg2rad(Dip)));...            % always downward on the fault plane (third component +)
           sin(deg2rad(Dip))];                                    % dip direction of the footwall

    Traction = Sigma_G * n_n;                                     % Traction acting on the fault in (N-E-d) system

    Sigma_n = transpose(n_n) * Traction;                          % normal stress on the fracture plane

    Tau_s = transpose(n_s) * Traction;                            % shear stress along strike
    Tau_d = transpose(n_d) * Traction;                            % shear stress along dip (+:downward, -:upward)

    Tau_vector = Tau_s*n_s + Tau_d*n_d;                           % total shear stress vector (N-E-down)
    Tau = norm(Tau_vector);                                       % total shear stress magnitude
end
