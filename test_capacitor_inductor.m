%example for how to use compute_planetary_motion(...)
function test_planetary_motion()
    C=10^-3;
    L=4.7*10^-4;
    i0 = 1;
    didt0 = 0;
    
    V0 = [i0; didt0];
    t_range = linspace(0,4*2*pi*sqrt(L*C),100);
    tspan = [t_range(1),t_range(end)];
    h_ref = 0.2;
    wrapper = @(t,V) [V(2);-V(1)/(C*L)];

    figure; 
    names = {"midpoint","kutta3rd","nystrom5th"};
    for i=1:length(names)
        BT_struct = rk_method(names{i});
        [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(wrapper,tspan,V0,h_ref,BT_struct);
        plot(t_list(1:4:end),X_list(1:4:end,1),'--','LineWidth', 1,"DisplayName",names{i});
        hold on
    end
    names = {"dormandprince","fehlberg","bogacki"};
    for i=1:length(names)
        BT_struct = rk_method(names{i});
        p = length(BT_struct.C)-1;
        error_desired = 10^-8;
        [t_list,X_list,h_avg, num_evals, failure_rate] = explicit_RK_variable_step_integration(wrapper,tspan,V0,h_ref,BT_struct,p,error_desired);
        plot(t_list(1:4:end),X_list(1:4:end,1),'--','LineWidth', 1,"DisplayName",names{i});
        hold on
    end
    legend(); 
end