run_global_truncation({@explicit_RK_fixed_step_integration, @explicit_RK_variable_step_integration},{"dormandprince"})

function run_global_truncation(step_funcs,methods)
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
    wrapper = @(t,V) gravity_rate_func(t,V,orbit_params);

    h_range = logspace(-3,0,100);
    e_range = logspace(-14,-8,100);
    iter = length(step_funcs)*length(methods);
    h_avg_list = zeros(iter,length(h_range));
    failure_list = zeros(iter,length(h_range));
    numerical = zeros(iter,length(h_range));
    evals = zeros(iter,length(h_range));
    for j = 1:iter
        meth = ceil(j/length(step_funcs));
        step = ceil(j/length(methods));
        step_func = step_funcs{step};
        BT_struct = rk_method(methods{meth});
        p = length(BT_struct.C)-1;
        if strcmp(func2str(step_func),func2str(@explicit_RK_fixed_step_integration))
            for i=1:length(h_range)
                [t_list,X_list,h_avg, num_evals] = step_func(wrapper,tspan,V0,h_range(i),BT_struct);
                h_avg_list(j,i) = h_avg;
                numerical(j,i) = norm(X_list(end, :)-V_list(end,:));
                evals(j,i) = num_evals;
            end
        end
        if strcmp(func2str(step_func),func2str(@explicit_RK_variable_step_integration))
            for i=1:length(e_range)
                [t_list,X_list,h_avg, num_evals,failure_rate] = step_func(wrapper,tspan,V0,h_range(1),BT_struct,p,e_range(i));
                h_avg_list(j,i) = h_avg;
                numerical(j,i) = norm(X_list(end, :)-V_list(end,:));
                evals(j,i) = num_evals;
                failure_list(j,i) = failure_rate;
            end
        end
    end
    figure;
    for j=1:iter
        meth = ceil(j/length(step_funcs));
        step = ceil(j/length(methods));
        filter_params = struct();
        filter_params.min_xval = 10^-2.5;
        filter_params.max_xval = 10^-1.5;
        loglog(h_avg_list(j,:),numerical(j,:),".",DisplayName=replace(func2str(step_funcs{step}),"_"," ")+" : "+methods{meth});
        hold on
        [p,k] = loglog_fit(h_avg_list(j,:),numerical(j,:),filter_params);
        loglog(h_avg_list(j,:),k*h_avg_list(j,:).^p,"-",DisplayName="p = "+string(p)+"; k = "+string(k));
    end
    xlabel("Step Size")
    ylabel("Global Truncation Error")
    axis([10^-3, 10^0,10^-12, 10^2]);
    legend("Location","northwest")  
    figure;
    for j=1:iter
        meth = ceil(j/length(step_funcs));
        step = ceil(j/length(methods));
        filter_params = struct();
        filter_params.min_xval = 10^3.5;
        loglog(evals(j,:),numerical(j,:),".",DisplayName=replace(func2str(step_funcs{step}),"_"," ")+" : "+methods{meth});
        hold on
        [p,k] = loglog_fit(evals(j,:),numerical(j,:),filter_params);
        loglog(evals(j,:),k*evals(j,:).^p,"-",DisplayName="p = "+string(p)+"; k = "+string(k));
    end
    xlabel("# of Calls")
    ylabel("Global Truncation Error")
    axis([10^2, 10^5, 10^-10, 10^2]);
    legend("Location","northeast")
    figure;
    for j=1:iter
        meth = ceil(j/length(step_funcs));
        step = ceil(j/length(methods));
        semilogx(h_avg_list(j,:),failure_list(j,:),".",DisplayName=replace(func2str(step_funcs{step}),"_"," ")+" : "+methods{meth});
        hold on
    end
    xlabel("Step Size")
    ylabel("Failure Rate")
    axis([10^-3, 10^0, 0, 0.1]);
    legend("Location","northeast")
end