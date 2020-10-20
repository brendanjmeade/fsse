function fault_system_parameter_sweep()
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
    dt = 0.0001;
    t = [0:dt:10];

    epsilon1_vec = linspace(0, 0.4, 100);
    epsilon2 = 0.1;
    ps_angle = 0.1 * cos(2 * t);
    dmat = zeros(numel(t), numel(epsilon1_vec));
    d1mat = zeros(numel(t), numel(epsilon1_vec));
    d2mat = zeros(numel(t), numel(epsilon1_vec));
    vvmat = zeros(numel(t), numel(epsilon1_vec));
    v1mat = zeros(numel(t), numel(epsilon1_vec));
    v2mat = zeros(numel(t), numel(epsilon1_vec));
    
    for i = 1:length(epsilon1_vec)
        epsilon1 = epsilon1_vec(i);
        [d, d1, d2, vv, v1, v2] = euler_integrate(t, v, [0, 0, 0], epsilon1, epsilon2, ...
                                      tau0, ps_angle, mu, H);
        dmat(:, i) = d;
        d1mat(:, i) = d1;
        d2mat(:, i) = d2;
        vvmat(:, i) = vv;
        v1mat(:, i) = v1;
        v2mat(:, i) = v2;
    end
    
    figure;    
    subplot(1, 2, 1)
    imagesc(v1mat)
    colormap(parula(3))
    xlabel("epsilon")
    ylabel("time")
    title("v1")
    colorbar;

    subplot(1, 2, 2)
    imagesc(v2mat)
    colormap(parula(3))
    xlabel("epsilon")
    ylabel("time")
    title("v2")
    colorbar;
    
end


function [d, d1, d2, vv, v1, v2] = euler_integrate(t, v, ics, epsilon1, epsilon2, ...
                                       tau0, ps_angle, mu, H)
    dt = t(2) - t(1);
    d = v * t;
    d1 = zeros(size(t));
    d2 = zeros(size(t));
    vv = v * ones(size(t));
    v1 = zeros(size(t));
    v2 = zeros(size(t));
 
    
    for i=2:length(t)
        sigma_xy = mu / H * (d(i-1) - d1(i-1) - d2(i-1));
        sigma_shear1(i) = sigma_xy * cos(ps_angle(i-1) + epsilon1);
        sigma_shear2(i) = sigma_xy * cos(ps_angle(i-1) + epsilon2);
        if sigma_shear1(i) >= tau0 && sigma_shear2(i) < tau0
            %increment fault1
            d1(i) = d1(i-1) + v * dt;
            d2(i) = d2(i-1);
            v1(i) = v;
            v2(i) = 0;
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) >= tau0
            %increment fault2
            d1(i) = d1(i-1);
            d2(i) = d2(i-1) + v * dt;
            v1(i) = 0;
            v2(i) = v;
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) < tau0
            %neither fault incremented
            d1(i) = d1(i-1);
            d2(i) = d2(i-1);
            v1(i) = 0;
            v2(i) = 0;
        else
            %split the increment with a 50% rule
            d1(i) = d1(i-1) + 0.5 * v * dt;
            d2(i) = d2(i-1) + 0.5 * v * dt;
            v1(i) = 0.5 * v;
            v2(i) = 0.5 * v;
        end
    end
end

