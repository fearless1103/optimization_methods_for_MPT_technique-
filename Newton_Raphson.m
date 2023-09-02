% Function to optimise the variables of dipole according to the receiving signals of magnetic field sensors 
for iter = 1:max_iteration                                                     %% number of iteration
    % Variables of dipole
    S = [sx, sy, sz];                                               %% Dipole position
    [ms,ns] = size(S);
    e = [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)];     %% Dipole direction
    vec_ep = e/norm(e);
    
    r_ps = P - repmat(S,mp/ms,np/ns);                               %% Distance between the dipole and sensors
    norm_r_ps = vecnorm(r_ps,2,2); 
    dot_e_r = dot(repmat(vec_ep,mp/ms,np/ns),r_ps,2);
    H = um/(4*pi)*(-vec_ep./norm_r_ps.^3 + 3*dot_e_r./norm_r_ps.^5.*r_ps);
    Hs = dot(H,vec_es,2);                                           %% Receiving signals of sensors
    
    % Approximate the derivatives by difference equation
    % Small increase of variables 
    sx_inc = sx*1.0001;
    sy_inc = sy*1.0001;
    sz_inc = sz*1.0001;
    phi_inc = phi*1.0001;
    theta_inc = theta*1.0001;
    % Corresponding receiving signals 
    Hs_inc_x = Hs_cal(sx_inc,sy,sz,phi,theta,P,vec_es);
    Hs_inc_y = Hs_cal(sx,sy_inc,sz,phi,theta,P,vec_es);
    Hs_inc_z = Hs_cal(sx,sy,sz_inc,phi,theta,P,vec_es);
    Hs_inc_phi = Hs_cal(sx,sy,sz,phi_inc,theta,P,vec_es);
    Hs_inc_theta = Hs_cal(sx,sy,sz,phi,theta_inc,P,vec_es);
    % Difference equation
    dHs_dx = (Hs_inc_x-Hs)/(sx*1e-4+1e-8);
    dHs_dy = (Hs_inc_y-Hs)/(sy*1e-4+1e-8);
    dHs_dz = (Hs_inc_z-Hs)/(sz*1e-4+1e-8);
    dHs_dphi = (Hs_inc_phi-Hs)/(phi*1e-4+1e-8);
    dHs_dtheta = (Hs_inc_theta-Hs)/(theta*1e-4+1e-8);
    
    % Variable update
    Lsdv = [dHs_dx, dHs_dy, dHs_dz, dHs_dphi, dHs_dtheta];                    % Jacobine/sensitivity matrix  
    S = Lsdv;
    rev = [sx; sy; sz; phi; theta];                                           % Variables of the current step       
    
    Hes = S' * S;                                                             % Approximate Hessian matrix
    [u1, s1, v1] = svd(Hes);
    s1 = diag(s1);
    s1 = diag(1./s1);   
    invH = v1*s1'*u1';                                                        % Inverse of the Hessian matrix
    
    rev = rev - invH * S' * (Hs - Hsm);                                       % Updated by Newton-Raphson
    
    % Limits the updated variables
    sx = min(max(rev(1),min(sx_a)),max(sx_a));
    sy = min(max(rev(2),min(sy_a)),max(sy_a));
    sz = min(max(rev(3),min(sz_a)),max(sz_a));
    phi = min(max(rev(4),min(phi_a)),max(phi_a));
    theta = min(max(rev(5),min(theta_a)),max(theta_a));
    
    res(iter)=sum(sum((abs(Hs-Hsm)).^2));
    
    Accuracy(iter,Ttest) = res(iter);
    Iteration(Ttest) = iter;
    
    % Threshold to stop according to residue
    if res(iter)<1e-6
        break
    end
end