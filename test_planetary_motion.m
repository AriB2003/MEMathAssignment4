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
    BT_struct.A = [0, 0; 0.5, 0];
    BT_struct.B = [0; 1];
    BT_struct.C = [0; 0.5];

    wrapper = @(t,V) gravity_rate_func(t,V,orbit_params);
    [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(wrapper,tspan,V0,h_ref,BT_struct);
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(X_list(:,1),X_list(:,2),'k');
    plot(V_list(:,1),V_list(:,2),'b');
end