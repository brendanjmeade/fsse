function [] = fault_system_idealized()
% Calculates displacements along 2 nearly parallel fault strands
% given some idealized assumptions:
% stress ~ sigma_xy + epsilon(t)
% fault1 ~ yield stress tau0, angle ~ 0 + epsilon1
% fault2 ~ yield stress tau0, angle ~ 0 + epsilon2
% Instead of changes directly to sigma_xy, below implements changes to
% principal stress angle as the perturbation
% yield criterion sigma_xy >= tau0
% System driven by applied constant velocity
% Stress is calculated elastically
% Written by VCT

sigma_xy = 1;
tau0 = 1;
v = 1;
mu = 1;
H = 1;
dt = 0.001;
t = [0:dt:10];
d = v*t;
d1 = 0*t;
d2 = d1;
sigma_shear1 = d1;
sigma_shear2 = d1;
% Redundant, but need initial d1,d2,sigma_shear1,sigma_shear2 all zero
d1(1) = 0;
d2(1) = 0;
sigma_shear1(1)=0;
sigma_shear2(1)=0;

% To toggle epsilon1 angle or ps_angle, pick different flags
eps1_flag = 2;
ps_angle_flag = 1;

if eps1_flag == 1
    epsilon1 = .05;
elseif eps1_flag == 2
    epsilon1 = .02;
else
    epsilon1 = 0;
end
epsilon2 = 0;
if ps_angle_flag == 1
    ps_angle = 0.1*cos(2*t); 
else
    ps_angle = 0*t;
end

for i=2:length(t)
    sigma_xy = mu/H*(d(i-1)-d1(i-1)-d2(i-1));
    sigma_shear1(i) = sigma_xy*cos(ps_angle(i-1)+epsilon1);
    sigma_shear2(i) = sigma_xy*cos(ps_angle(i-1)+epsilon2);
    if sigma_shear1(i) >= tau0 && sigma_shear2(i) < tau0
        %increment fault1
        d1(i)=d1(i-1)+v*dt;
        d2(i)=d2(i-1);
    elseif sigma_shear1(i) < tau0 && sigma_shear2(i) >= tau0
        %increment fault2
        d1(i)=d1(i-1);
        d2(i)=d2(i-1)+v*dt;
    elseif sigma_shear1(i) < tau0 && sigma_shear2(i) < tau0
        %neither fault incremented
        d1(i)=d1(i-1);
        d2(i)=d2(i-1);
    else
        %split the increment with some rule
        % 50% rule
        d1(i)=d1(i-1)+0.5*v*dt;
        d2(i)=d2(i-1)+0.5*v*dt;
    end
end

plot(t,d,t,d1,t,d2,t,sigma_shear1,t,sigma_shear2);
ylim([-1 10]);
xlabel('time');
legend('d total','d1','d2','stress 1','stress 2');