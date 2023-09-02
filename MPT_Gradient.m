clear all; close all; clc

%% Definition
% Sensor positions (16 sensors arranged in 4 layers )
% p in the paper
pr = 10e-3;                                                         %% Radius of sensors
p_th_int = pi/2;                                                    %% Central angle interval of sensors
p_x = repmat(pr*cos(0:p_th_int:2*pi-p_th_int)',4,1);                %% x positions of sensors
p_y = repmat(pr*sin(0:p_th_int:2*pi-p_th_int)',4,1);                %% y positions of sensors
p_z = [0*ones(1,4), 7*ones(1,4) 14*ones(1,4) 21*ones(1,4)]'*1e-3;   %% z positions of sensors
P = [p_x, p_y, p_z];                                                %% Coordinates of sensors center
[mp,np] = size(P);

% Sensor normal vector
% e_s in the paper
beta = angle(p_x+1j*p_y);                                           %% Horizontal polar angle of center of sensors
alpha = angle(pr*ones(16,1)+1j*p_z);                                %% Longitudinal polar angle of center of sensors
es = [cos(beta), sin(beta), sin(alpha)];                            %% Normal vector of sensors
vec_es = es./vecnorm(es,2,2);

%% Limits of possible values of the dipole
sx_a = [-pr pr];                                                   %% x position
sy_a = sx_a;                                                        %% y position
sz_a = [min(p_z) max(p_z)];                                         %% z position
phi_a = [-2*pi 2*pi];                                               %% Horizontal polar angle of the dipole
theta_a = [-pi pi];                                                 %% Longitudinal polar angle of the dipole

%% Simulate the magnetic field signals of sensors
%actual dipole position for Gradient Descent method
sx_real = 2e-3;
sy_real = 3e-3;
sz_real = 3e-3;
S = [sx_real, sy_real, sz_real];
[ms,ns] = size(S);
phi_real = pi/2;
theta_real = pi/3;

% Normal vector of dipole pivot (e_p in the paper)
e = [cos(phi_real)*cos(theta_real), sin(phi_real)*cos(theta_real), sin(theta_real)];
vec_ep = e/norm(e);
% Distance vector between the center of dipole and sensors (r_ps in the paper)
r_ps = P - repmat(S,mp/ms,np/ns);
um = 4*pi*1e-7;

%% H Calculation
norm_r_ps = vecnorm(r_ps,2,2);                                      %% Norm of distance vector
dot_e_r = dot(repmat(vec_ep,mp/ms,np/ns),r_ps,2);                   %% dot product of e_p .* r_ps
H = um/4/pi*(-vec_ep./norm_r_ps.^3 + 3*dot_e_r./norm_r_ps.^5.*r_ps);%% H field at the center of sensors
Hsm = dot(H,vec_es,2);                                              %% Measured H field along the sensitive axis of sensors

%% Inverse
max_test = 10;
max_iteration = 50;
Accuracy = zeros(max_iteration,max_test);
dipole_record = [];
for Ttest = 1:max_test
  % random initial values for gradient descent method  
    sx = unifrnd(-pr,pr);
    sy = unifrnd(-pr,pr);
    sz = unifrnd(min(p_z),max(p_z));
    phi = 2*pi/randi(9);
    theta = 2*pi/randi(9);
    dipole_record(:,Ttest) = [sx,sy,sz,phi,theta];
    
    % Optimization based on the Gradient Descent method
    % Step size (learning rate) for Gradient Descent
    alpha = 1e-5;

    tStart = tic;
    for iter = 1:max_iteration
        S = [sx, sy, sz];
        [ms,ns] = size(S);
        e = [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)];
        vec_ep = e / norm(e);
        
        r_ps = P - repmat(S, mp/ms, np/ns);
        norm_r_ps = vecnorm(r_ps, 2, 2);
        dot_e_r = dot(repmat(vec_ep, mp/ms, np/ns), r_ps, 2);
        H = um / (4 * pi) * (-vec_ep ./ norm_r_ps.^3 + 3 * dot_e_r ./ norm_r_ps.^5 .* r_ps);
        Hs = dot(H, vec_es, 2);
        
        residue = Hs - Hsm;
        
        gradient = zeros(1, 5);
        for dim = 1:5
            % Small perturbation
            perturbation = zeros(1, 5);
            perturbation(dim) = 0.0001;
            
            Hs_perturbed = Hs_cal(sx + perturbation(1), sy + perturbation(2), ...
                sz + perturbation(3), phi + perturbation(4), ...
                theta + perturbation(5), P, vec_es);
            gradient(dim) = sum((Hs_perturbed - Hs) .* perturbation(dim)) ./ ((perturbation(dim) + 1e-8));
        end
        
        % Update dipole variables
        sx = sx - alpha * sum(residue .* (-gradient(1))) / (norm(gradient(1))^2 + 1e-8);
        sy = sy - alpha * sum(residue .* (-gradient(2))) / (norm(gradient(2))^2 + 1e-8);
        sz = sz - alpha * sum(residue .* (-gradient(3))) / (norm(gradient(3))^2 + 1e-8);
        phi = phi - alpha * sum(residue .* (-gradient(4))) / (norm(gradient(4))^2 + 1e-8);
        theta = theta - alpha * sum(residue .* (-gradient(5))) / (norm(gradient(5))^2 + 1e-8);
        
        res(iter) = sum(sum((abs(Hs - Hsm)).^2));
        
        Accuracy(iter,Ttest) = res(iter);
        Iteration(Ttest) = iter;
        
        % Threshold to stop according to residue
        if res(iter) < 1e-6
            break
        end
    end
    dipole_result(:,Ttest) = [sx,sy,sz,phi,theta];
    Time(Ttest) = toc(tStart);
    
end
