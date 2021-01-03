function fault_system_epsilon2_viscous_dimensions()
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
    close all;
    tau0 = 1;
    v = 1;
    mu = 1;
    H = 1;
    dt = 0.01;
    t = [0:dt:10];
        
    epsilon1 = 0.05;
    epsilon2 = 0.0;
    [d_total_plastic, d_1_plastic, d_2_plastic] = euler_integrate(t, v, [0, 0, 0], epsilon1, epsilon2, tau0, mu, H);
    [d_total_viscous, d_1_viscous, d_2_viscous] = euler_integrate_viscous(t, v, [0, 0, 0], epsilon1, epsilon2, tau0, mu, H);

    figure;
    hold on;
    plot(t, d_total_plastic, "-k")
    plot(t, d_1_plastic, "-r")
    plot(t, d_2_plastic, "-b")
    % plot(t, d_total_viscous, "--k")
    % plot(t, d_1_viscous, "--r")
    % plot(t, d_2_viscous, "--b")
    xlabel("time (-)");
    box on;
    legend("d", "d1", "d2");
end


function [d, d1, d2] = euler_integrate(t, v, ics, epsilon1, epsilon2, tau0, mu, H)
    dt = t(2) - t(1);
    d = v * t;
    d1 = zeros(size(t));
    d2 = zeros(size(t));
    
    stress_scale_factor = 2e0; % Pa - FOR DIMENSIONALIZATION
    tau0 = stress_scale_factor * tau0;
    for i=2:length(t)
        sigma_xy = mu / H * (d(i-1) - d1(i-1) - d2(i-1));
        sigma_xy = stress_scale_factor * sigma_xy;
        
        % Victor's new lines from late October
        sigma_shear1(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon1)+sigma_xy*cos(2*epsilon1);
        sigma_shear2(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon2)+sigma_xy*cos(2*epsilon2);        
        disp([sigma_shear1(i) / tau0, sigma_shear2(i) / tau0])
        
        if sigma_shear1(i) >= tau0 && sigma_shear2(i) < tau0
            % increment fault 1
            d1(i) = d1(i-1) + v * dt;
            d2(i) = d2(i-1);
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) >= tau0
            % increment fault 2
            d1(i) = d1(i-1);
            d2(i) = d2(i-1) + v * dt;
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) < tau0
            % neither fault incremented
            d1(i) = d1(i-1);
            d2(i) = d2(i-1);
        else
            % split the increment with a 50% rule
            d1(i) = d1(i-1) + 0.5 * v * dt;
            d2(i) = d2(i-1) + 0.5 * v * dt;
        end
    end
end


function [d, d1, d2] = euler_integrate_viscous(t, v, ics, epsilon1, epsilon2, tau0, mu, H)
    dt = t(2) - t(1);
    d = v * t;
    d1 = zeros(size(t));
    d2 = zeros(size(t));
            
    for i=2:length(t)
        sigma_xy = mu / H * (d(i-1) - d1(i-1) - d2(i-1));

        % Victor's new lines from late October
        sigma_shear1(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon1)+sigma_xy*cos(2*epsilon1);
        sigma_shear2(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon2)+sigma_xy*cos(2*epsilon2);

        length_scale = 1;
        eta = 1;
        n = 0.001;
        v1 = sign(sigma_shear1(i))*(abs(sigma_shear1(i))^(1/n)) * length_scale / eta;
        v2 = sign(sigma_shear1(i))*(abs(sigma_shear2(i))^(1/n)) * length_scale / eta;
        vtotal = v1 + v2;
        % v1 = v1 / vtotal;
        % v2 = v2 / vtotal;
        d1(i) = d1(i-1) + v1 * dt;
        d2(i) = d2(i-1) + v2 * dt;
    end
end
