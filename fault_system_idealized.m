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

    % Store fixed, specified time, and non-updating parameters in global structure
    % for easy access in ODE solver
    global G
    G.tau0 = tau0;
    G.mu = mu;
    G.H = H;
    G.v = v;
    G.epsilon1 = epsilon1;
    G.epsilon2 = epsilon2;
    G.ps_angle_flag = ps_angle_flag;
    
    % ODE initial conditions and ode45 solve
    ics = [0; 0; 0];
    tspan = [min(t) max(t)];
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t45, y45] = ode45(@derivs, tspan, ics, opts);
    % Implicit stiff solver does not get past t = 1
    % [t15s, y15s] = ode15s(@derivs, tspan, ics, opts);
    
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

    figure;
    hold on;
    plot(t45, y45(:, 1), "-k");
    plot(t45, y45(:, 2), "-r");
    plot(t45, y45(:, 3), "-b");
    plot(t, d, "--k")
    plot(t, d1, "--r")
    plot(t, d2, "--b")
    xlabel('time');
    box on;
    legend('d total - ode45', 'd1 - ode45', 'd2 - ode45', 'd total - Euler', 'd1 - Euler', 'd2 - Euler');
    
    figure;
    start_idx = 1000000;
    subplot(2, 2, 1);
    plot(d1(start_idx:end), sigma_shear1(start_idx:end), 'r.');
    xlabel("displacment");
    ylabel("shear stress");
    title("displacement fault 1 - shear stress fault 1");
    
    subplot(2, 2, 2);
    plot(d1(start_idx:end), sigma_shear2(start_idx:end), 'r.');
    xlabel("displacment");
    ylabel("shear stress");
    title("displacement fault 1 - shear stress fault 2");

    subplot(2, 2, 3);
    plot(d2(start_idx:end), sigma_shear1(start_idx:end), 'r.');
    xlabel("displacment");
    ylabel("shear stress");
    title("displacement fault 2 - shear stress fault 1");

    subplot(2, 2, 4);
    plot(d2(start_idx:end), sigma_shear2(start_idx:end), 'r.');
    xlabel("displacment");
    ylabel("shear stress");
    title("displacement fault 2 - shear stress fault 2");
    
end


function dy_dt = derivs(t, y)
    global G;    
    d = y(1);
    d1 = y(2);
    d2 = y(3);

    if G.ps_angle_flag == 1
        ps_angle = 0.1*cos(2*t); 
    else
        ps_angle = 0*t;
    end

    sigma_xy = G.mu / G.H * (d - d1 - d2);
    sigma_shear1 = sigma_xy * cos(ps_angle + G.epsilon1);
    sigma_shear2 = sigma_xy * cos(ps_angle + G.epsilon2);

    % Derivatives
    d_dt = G.v;
    if sigma_shear1 >= G.tau0 && sigma_shear2 < G.tau0
        % Increment fault 1
        d1_dt = G.v;
        d2_dt = 0;
    elseif sigma_shear1 < G.tau0 && sigma_shear2 >= G.tau0
        % Increment fault 2
        d1_dt = 0;
        d2_dt = G.v;
    elseif sigma_shear1 < G.tau0 && sigma_shear2 < G.tau0
        % Neither fault incremented
        d1_dt = 0;
        d2_dt = 0;
    else
        % Split the increment with some rule
        % 50% rule
        d1_dt = 0.5 * G.v;
        d2_dt = 0.5 * G.v;
    end
    dy_dt = [d_dt; d1_dt; d2_dt];    
end

