run_local_truncation({@explicit_RK_step},{"midpoint","kutta3rd","nystrom5th"})

function run_local_truncation(step_funcs, methods)
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    wrapper = @(t,V) gravity_rate_func(t,V,orbit_params);

    h_range = logspace(-5,1,100);
    iter = length(step_funcs)*length(methods);
    analytical = zeros(size(h_range,1));
    for i=1:length(h_range)
        analytical(i) = norm(compute_planetary_motion(0+h_range(i),V0,orbit_params)-V0);
    end
    numerical = zeros(iter,length(h_range));
    for j = 1:iter
        for i=1:length(h_range)
            meth = ceil(j/length(step_funcs));
            step = ceil(j/length(methods));
            XB = compute_planetary_motion(0+h_range(i),V0,orbit_params);
            step_func = step_funcs{step};
            [predicted_XB,num_evals] = step_func(wrapper,0,V0,h_range(i),rk_method(methods{meth}));
            numerical(j,i) = norm(predicted_XB-XB);
        end
    end
    figure;
    loglog(h_range,analytical,".",DisplayName="analytical");
    hold on
    filter_params = struct();
    filter_params.min_xval = 10^-2;
    filter_params.max_xval = 10^0;
    [p,k] = loglog_fit(h_range,analytical,filter_params);
    loglog(h_range,k*h_range.^p,"-",DisplayName="p = "+string(p)+"; k = "+string(k));
    for j=1:iter
        meth = ceil(j/length(step_funcs));
        step = ceil(j/length(methods));
        loglog(h_range,numerical(j,:),".",DisplayName=replace(func2str(step_funcs{step}),"_"," ")+" : "+methods{meth});
        hold on
        [p,k] = loglog_fit(h_range,numerical(j,:),filter_params);
        loglog(h_range,k*h_range.^p,"-",DisplayName="p = "+string(p)+"; k = "+string(k));
    end
    xlabel("Step Size")
    ylabel("Local Truncation Error")
    axis([10^-5, 10^1,10^-16, 10^2]);
    legend("Location","northwest")
end