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
    h_ref = 0.1;
    BT_struct = struct();
    BT_struct.A = [0, 0, 0, 0, 0, 0; 
                   1/4, 0, 0, 0, 0, 0;
                   4/25, 6/25, 0, 0, 0, 0;
                   1/4, -3, 15/4, 0, 0, 0;
                   2/27, 10/9, -50/81, 8/81, 0, 0;
                   2/25, 12/25, 2/15, 8/75, 0, 0];
    BT_struct.B = [23/192; 0; 125/192; 0; -27/64; 125/192];
    BT_struct.C = [0; 1/3; 2/5; 1; 2/3; 4/5];

    wrapper = @(t,V) gravity_rate_func(t,V,orbit_params);
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(wrapper,tspan,V0,h_ref,BT_struct);
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'b','linewidth',2);
    plot(X_list(1:4:end,1),X_list(1:4:end,2),'ko','markersize', 3);
    

end