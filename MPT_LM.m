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
sx_a = [-pr pr];                                                    %% x position
sy_a = sx_a;                                                        %% y position
sz_a = [min(p_z) max(p_z)];                                         %% z position
phi_a = [-2*pi 2*pi];                                               %% Horizontal polar angle of the dipole
theta_a = [-pi pi];                                                 %% Longitudinal polar angle of the dipole

%% Simulate the magnetic field signals of sensors
% give a certain dipole position
sx_real = 2e-3;                 
sy_real = 3e-3;
sz_real = 3e-3;
S = [sx_real, sy_real, sz_real];
[ms,ns] = size(S);
phi_real = pi/2;
theta_real = pi/3;

% Normal vector of dipole pivot (e_p in the paper)
e = [cos(phi_real) * cos(theta_real), sin(phi_real) * cos(theta_real), sin(theta_real)];
vec_ep = e / norm(e);
% Distance vector between the center of dipole and sensors (r_ps in the paper)
r_ps = P - repmat(S, mp / ms, np / ns);
um = 4 * pi * 1e-7;

%% H Calculation
norm_r_ps = vecnorm(r_ps, 2, 2);
dot_e_r = dot(repmat(vec_ep, mp / ms, np / ns), r_ps, 2);
H = um / (4 * pi) * (-vec_ep ./ norm_r_ps.^3 + 3 * dot_e_r ./ norm_r_ps.^5 .* r_ps);
Hs = dot(H, vec_es, 2);

%% Levenberg-Marquardt Optimization
max_test = 10;
max_iteration = 20;
Accuracy = zeros(max_iteration, max_test);
dipole_record = [];

for Ttest = 1:max_test
    %initial values randomly produced
    sx = unifrnd(-pr,pr);
    sy = unifrnd(-pr,pr);
    sz = unifrnd(min(p_z),max(p_z));
    phi = 2*pi/randi(9);
    theta = 2*pi/randi(9);
    dipole_record(:,Ttest) = [sx,sy,sz,phi,theta];
    
    resnorm = [];
    tStart = tic;
    for iter = 1:max_iteration
        options = optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12);
        [optimized_params, ~, residual] = lsqnonlin(@(params) calc_residual(params, P, vec_es, Hs), ...
            [sx, sy, sz, phi, theta], [], [], options);
        
        %Calculate the Jacobian matrix
        damping_factor = 1e-5;
        J = calc_jacobian(optimized_params, P, vec_es, Hs);
        
        delta_params = (J' * J + damping_factor * diag(diag(J' * J))) \ (J' * residual);
        
        %Update parameters
        optimized_params = optimized_params + delta_params;
        sx = optimized_params(1);
        sy = optimized_params(2);
        sz = optimized_params(3);
        phi = optimized_params(4);
        theta = optimized_params(5);
        
        resnorm = norm(residual,2);
        Accuracy(iter,Ttest) = resnorm;
        
        if resnorm < 1e-6
            break
        end
    end
    
    sx_final = optimized_params(1,1);
    sy_final = optimized_params(1,2);
    sz_final = optimized_params(1,3);
    phi_final = optimized_params(1,4);
    theta_final = optimized_params(1,5);

    Time(Ttest) = toc(tStart);
    dipole_result(:,Ttest) = [sx_final,sy_final,sz_final,phi_final,theta_final];
end

function F = calc_residual(params, P, vec_es, Hs)
% Calculate the difference between measured and simulated signals
S = [params(1), params(2), params(3)];
e = [cos(params(4)) * cos(params(5)), sin(params(4)) * cos(params(5)), sin(params(5))];
vec_ep = e / norm(e);

r_ps = P - repmat(S, size(P, 1), 1);
um = 4 * pi * 1e-7;

norm_r_ps = vecnorm(r_ps, 2, 2);
dot_e_r = dot(repmat(vec_ep, size(P, 1), 1), r_ps, 2);

H = um / (4 * pi) * (-vec_ep ./ norm_r_ps.^3 + 3 * dot_e_r ./ norm_r_ps.^5 .* r_ps);
Hs_simulated = dot(H, vec_es, 2);
F = Hs - Hs_simulated;
end

function J = calc_jacobian(params, P, vec_es, Hs)
    % Calculate the Jacobian matrix
    delta = 1e-8;
    num_params = length(params);
    num_sensors = size(P, 1);
    J = zeros(num_sensors, num_params);

    for i = 1:num_params
        perturbed_params1 = params;
        perturbed_params2 = params;
        perturbed_params1(i) = perturbed_params1(i) - delta/2;
        perturbed_params2(i) = perturbed_params2(i) + delta/2;
        
        F_perturbed1 = calc_residual(perturbed_params1, P, vec_es, Hs);
        F_perturbed2 = calc_residual(perturbed_params2, P, vec_es, Hs);
        
        J(:, i) = (F_perturbed2 - F_perturbed1) / delta;
    end
end