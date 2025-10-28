%example for how to use compute_planetary_motion(...)
function test_planetary_motion()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

    tspan = [t_range(1),t_range(end)];
    h_ref = 0.2;
    wrapper = @(t,V) gravity_rate_func(t,V,orbit_params);

    figure; 
    axis equal; axis square;
    axis([-3,9,-6,6])
    hold on
    plot(0,0,'ro','markerfacecolor','r','MarkerSize',10, "DisplayName","sun");
    plot(V_list(:,1),V_list(:,2),'b', "DisplayName","analytical");

    names = {"midpoint","kutta3rd","nystrom5th"};
    for i=1:length(names)
        BT_struct = rk_method(names{i});
        [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(wrapper,tspan,V0,h_ref,BT_struct);
        plot(X_list(1:4:end,1),X_list(1:4:end,2),'--','LineWidth', 1,"DisplayName",names{i});
    end
    names = {"dormandprince","fehlberg","bogacki"};
    for i=1:length(names)
        BT_struct = rk_method(names{i});
        p = length(BT_struct.C)-1;
        error_desired = 10^-8;
        [t_list,X_list,h_avg, num_evals, failure_rate] = explicit_RK_variable_step_integration(wrapper,tspan,V0,h_ref,BT_struct,p,error_desired);
        plot(X_list(1:4:end,1),X_list(1:4:end,2),'--','LineWidth', 1,"DisplayName",names{i});
    end
    legend();  
    figure;
    inc = 25;
    subplot(2,1,1);
    plot(t_list(1:inc:end), X_list(1:inc:end,1),'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2, DisplayName="X")
    hold on
    plot(t_list(1:inc:end), X_list(1:inc:end,2),'go-','markerfacecolor','k','markeredgecolor','k','markersize',2, DisplayName="Y")
    ylabel("Position (m)")
    legend()
    subplot(2,1,2);
    plot(t_list(1:inc:end), X_list(1:inc:end,3),'ro-','markerfacecolor','k','markeredgecolor','k','markersize',2, DisplayName="Vx")
    hold on
    plot(t_list(1:inc:end), X_list(1:inc:end,4),'go-','markerfacecolor','k','markeredgecolor','k','markersize',2, DisplayName="Vy")
    ylabel("Velocity (m/s)")
    xlabel("Time (s)")
    legend()
    h_list = diff(t_list);
    r = sqrt(X_list(:,1).^2+X_list(:,2).^2);
    figure;
    semilogy(r(2:end),h_list,".")
    xlabel("Distance (m)")
    ylabel("Step Size (s)")
end