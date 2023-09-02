% Function to calculate the receiving signals of magnetic field sensors for a given dipole
function Hs = Hs_cal(sx,sy,sz,phi,theta,P,vec_es)
    [mp,np] = size(P);
    S = [sx, sy, sz];
    [ms,ns] = size(S);
    e = [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)];
    vec_ep = e/norm(e);
    um = 4*pi*1e-7;
    
    % Receiving signal calculation
    r_ps = P - repmat(S,mp/ms,np/ns);
    norm_r_ps = vecnorm(r_ps,2,2);
    dot_e_r = dot(repmat(vec_ep,mp/ms,np/ns),r_ps,2);
    H = um/4/pi*(-vec_ep./norm_r_ps.^3 + 3*dot_e_r./norm_r_ps.^5.*r_ps);
    Hs = dot(H,vec_es,2);

end
