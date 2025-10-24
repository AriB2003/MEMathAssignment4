%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)

 num_evals = 0;

    % calculate the timesteps at which the numerical diff. method will be run
    t_delta = tspan(2)-tspan(1);
    N = ceil(t_delta/h_ref);
    h_avg = t_delta/N;
    %t_list = (tspan(1):h_avg:tspan(2))';
    t_list = [tspan(1)];

    % initialize storage structure for the X values
    X_list = X0';
    
    h_next = h_avg;
    i = 2;

    % run numerical differentiation method based on parameters given above
    while t_list(i-1) < tspan(2)
        redo = true;
        t_list = [t_list,t_list(i-1)+h_next];
        while redo == true
          [X_next, evals, h_next, redo] = explicit_RK_variable_step(rate_func_in, t_list(i-1), X_list(i-1,:)', h_next, BT_struct, p , error_desired);  
          num_evals = num_evals + evals;
        end  
        X_list = [X_list;X_next'];
        i = i+1;
    end

end