clear all; close all; clc

%% Definition
% Sensor positions (16 sensors arranged in 4 layers )
% p in the paper
pr = 10e-3;                                                         %% Radius of sensors
p_th_int = pi/2;                                                    %% Central angle interval of sensors 
p_x = repmat(pr*cos(0:p_th_int:2*pi-p_th_int)',4,1);                            %% x positions of sensors
p_y = repmat(pr*sin(0:p_th_int:2*pi-p_th_int)',4,1);                            %% y positions of sensors
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
% give a certain dipole position
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
                                                                    % along the normal vector of sensors
%% Inverse
% % specific initial values of dipole 
% sx = 1e-3; 
% sy = 1e-3;
% sz = 1e-3;
% phi = pi/6;
% theta = pi/6;
% max_iteration = 20;
% run Newton_Raphson.m

max_iteration = 20;
max_test = 10;
Accuracy = zeros(max_iteration,max_test); 
dipole_record = [];
for Ttest = 1:max_test
    %initial values randomly produced%
    sx = unifrnd(-pr,pr);
    sy = unifrnd(-pr,pr);
    sz = unifrnd(min(p_z),max(p_z));
    phi = 2*pi/randi(9);
    theta = 2*pi/randi(9);
    dipole_record(:,Ttest) = [sx,sy,sz,phi,theta];
    
    %Processing time for optimization process%
    tStart = tic;
    
    run Newton_Raphson.m
    % Optimisation based on the Newton-Raphson method
    dipole_result(:,Ttest) = [sx,sy,sz,phi,theta];
    
    Time(Ttest) = toc(tStart);
end

%% inverse results 
% Calculate the relative error of variables
eps_re_x = abs((sx-sx_real)/(sx_real+1e-8))*100;
eps_re_y =  abs((sy-sy_real)/(sy_real+1e-8))*100;
eps_re_z =  abs((sz-sz_real)/(sz_real+1e-8))*100;
eps_re_phi =  abs((phi-phi_real)/(phi_real+1e-8))*100;
eps_re_theta =  abs((theta-theta_real)/(theta_real+1e-8))*100;

% Information output
disp(['The actual x postion is ' num2str(sx_real*1e3) '(mm) and the estimated value is ' num2str(sx*1e3) '(mm)'])
disp(['The actual y postion is ' num2str(sy_real*1e3) '(mm) and the estimated value is ' num2str(sy*1e3) '(mm)'])
disp(['The actual z postion is ' num2str(sz_real*1e3) '(mm) and the estimated value is ' num2str(sz*1e3) '(mm)'])
disp(['The actual angle \Phi is ' num2str(phi_real) '(rad/s) and the estimated value is ' num2str(phi) '(rad/s)'])
disp(['The actual angle \Theta is ' num2str(theta_real) '(rad/s) and the estimated value is ' num2str(theta) '(rad/s)'])
disp(['The relative error of the x position is  ' num2str(eps_re_x) '%'])
disp(['The relative error of the y position is  ' num2str(eps_re_y) '%'])
disp(['The relative error of the z position is  ' num2str(eps_re_z) '%'])
disp(['The relative error of the angle \Phi is  ' num2str(eps_re_phi) '%'])
disp(['The relative error of the angle \Theta is  ' num2str(eps_re_theta) '%'])

% Measured magnetic field on sensors
figure,plot(Hs,'*');                                                %% Given signals
hold on 
plot(Hsm,'o');                                                      %% Optimised signals
% set(gca,'xscale','log')
xlabel('Sensor index','interpreter','latex', 'Fontname','Times New Roman','fontsize',18)
ylabel('$H_s$ (T)','interpreter','latex', 'Fontname','Times New Roman','fontsize',18)
set(gcf,'position',[200 200 1100 400]);
legend('Optimised','Measurement')

% Optimisation residual
figure,plot(res)
xlabel('Iterative step','interpreter','latex', 'Fontname','Times New Roman','fontsize',18)
ylabel('Residual','interpreter','latex', 'Fontname','Times New Roman','fontsize',18)

%% show position
clear alpha
clear figure(2)
pic = figure(2);
res = 16;
vx = pr*[-1+2/res:2/res:1];
vy = pr*[-1+2/res:2/res:1];
vz = max(p_z)*[-1+2/res:2/res:1];
[x,y,z] = meshgrid(vx*1e3,vy*1e3,vz*1e3);

ys = vy*1e3;
xs = vx*1e3;
zs = vz*1e3;
vol_pos = zeros(length(vx),length(vy),length(vz));
pos_x = round((sx-min(sx_a))/(max(sx_a)-min(sx_a))*length(vx));
pos_y = round((sy-min(sy_a))/(max(sy_a)-min(sy_a))*length(vy));
pos_z = round((sz-min(sz_a))/(max(sz_a)-min(sz_a))*length(vz));
vol_pos(pos_x,pos_y,pos_z+1) = 1;
c = vol_pos;
h = slice(x,y,z,c,xs,ys,zs,'cubic');
set(h,'FaceColor','interp','EdgeColor','interp')
set(gcf, 'Color', 'w');
camproj perspective
box on
view(-200,200)
colormap(flipud(pink))
set(gca, 'CLim', [max(-0.1*max(c,[],"all"),0), min(max(c,[],"all")*0.1,100)]);
hcb = colorbar;
cbh = findall(pic, 'Type', 'ColorBar');
cTH = get(cbh,'Title');

ylabel('$y$ (mm)','interpreter','latex', 'Fontname','Times New Roman','fontsize',20)
zlabel('$z$ (mm)','interpreter','latex', 'Fontname','Times New Roman','fontsize',20)
xlabel('$x$ (mm) ','interpreter','latex', 'Fontname','Times New Roman','fontsize',20)
set(gcf,'position',[200 200 460 400]);
grid off 
shading interp
set(gcf, 'PaperPosition', [-2.05 -1.75 22 22]);
set(gcf, 'PaperSize', [20 20]); 
alpha(0.07) 
view([-1,-1,0.7])
set(hcb,'Position',[0.90 0.2 0.03 0.55]);
set(gca,'FontSize',20);