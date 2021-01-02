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
    % Original values
    % tau0 = 1;
    % v = 1;
    % mu = 1;
    % H = 1;
    % dt = 0.001;
    % t = [0:dt:10];
    
    % Dimensionalized values
    siay = 60 * 60 * 24 * 365.25; % seconds in a year
    tau0 = 15e6; % critical switching stress (Pa)
    v = 0.05 / siay; % reference velocity (m/s) - 50 mm/yr
    mu = 3e10; % shear modulus (Pa)
    H = 50e3; % length scale (km)
    t = [0:0.1:1e4]; % (years!)
    t = siay * t;
    
    epsilon1 = 0.05;
    epsilon2 = 0.0;
    ps_angle = 0.1 * cos(2 * t); 
    [d_total_plastic, d_1_plastic, d_2_plastic] = euler_integrate(t, v, [0, 0, 0], epsilon1, epsilon2, tau0, ps_angle, mu, H);
    [d_total_viscous, d_1_viscous, d_2_viscous] = euler_integrate_viscous(t, v, [0, 0, 0], epsilon1, epsilon2, tau0, ps_angle, mu, H);
    
    figure;
    hold on;
    plot(t/siay, d_total_plastic, "-k")
    plot(t/siay, d_1_plastic, "-r")
    plot(t/siay, d_2_plastic, "-b")
    % plot(t, d_total_viscous, "--k")
    % plot(t, d_1_viscous, "--r")
    % plot(t, d_2_viscous, "--b")
    xlabel("time (years)");
    box on;
    legend("d", "d1", "d2");
end


function [d, d1, d2] = euler_integrate(t, v, ics, epsilon1, epsilon2, tau0, ps_angle, mu, H)
    dt = t(2) - t(1);
    d = v * t;
    d1 = zeros(size(t));
    d2 = zeros(size(t));
            
    for i=2:length(t)
        sigma_xy = mu / H * (d(i-1) - d1(i-1) - d2(i-1));

        % Victor's new lines from late October
        % sigma_shear1(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon1)+sigma_xy*cos(2*epsilon1);
        % sigma_shear2(i) = 0.1*cos(2*t(i-1))*sin(2*epsilon2)+sigma_xy*cos(2*epsilon2);
        % sigma_shear1(i) = 6e3*cos(2*t(i-1))*sin(2*epsilon1)+sigma_xy*cos(2*epsilon1);
        % sigma_shear2(i) = 6e3*cos(2*t(i-1))*sin(2*epsilon2)+sigma_xy*cos(2*epsilon2);
        sigma_shear1(i) = 2e5 * cos(2 * t(i-1)) * sin(2 * epsilon1) + sigma_xy * cos(2 * epsilon1);
        sigma_shear2(i) = 2e5 * cos(2 * t(i-1)) * sin(2 * epsilon2) + sigma_xy * cos(2 * epsilon1);
        
        % [sigma_shear1(i) / tau0, sigma_shear2(i) / tau0]
        
        if sigma_shear1(i) >= tau0 && sigma_shear2(i) < tau0
            % increment fault 1
            disp("Case 1")
            d1(i) = d1(i-1) + v * dt;
            d2(i) = d2(i-1);
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) >= tau0
            % increment fault 2
            disp("Case 2")
            d1(i) = d1(i-1);
            d2(i) = d2(i-1) + v * dt;
        elseif sigma_shear1(i) < tau0 && sigma_shear2(i) < tau0
            % neither fault incremented
            disp("Case 3")
            d1(i) = d1(i-1);
            d2(i) = d2(i-1);
        else
            % split the increment with a 50% rule
            disp("Case 4")
            d1(i) = d1(i-1) + 0.5 * v * dt;
            d2(i) = d2(i-1) + 0.5 * v * dt;
        end
    end
end


function [d, d1, d2] = euler_integrate_viscous(t, v, ics, epsilon1, epsilon2, tau0, ps_angle, mu, H)
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
