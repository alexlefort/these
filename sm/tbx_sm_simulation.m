function res = tbx_sm_simulation(param, noise , mission, pctrl, is_saved)
    
    rng(1234);
    close all;
    
    %% Initial State
    
    state  = mission.state_init ;
    cmd    = mission.cmd_init   ;
    intern = pctrl.intern_init  ;
    res    = tbx_sm_make_simu_vect(mission.T);
    
    %% Simulation Loop
    t = 0.0;
    
    for ii=1:length(mission.T)
         t = t + pctrl.dt;
         state        = tbx_sm_dynamics(state, cmd, param, pctrl.dt);
         measure      = tbx_sm_sensors(state, noise);
         
         order.ZCo    = mission.ZCo(ii);
         order.PsiCo  = mission.PsiCo(ii);
         order.Prop   = mission.Prop(ii);
         
         [cmd,intern] = tbx_sm_controller(measure, order, pctrl, intern);
         res          = tbx_sm_save_simu_vect(res, state, measure, cmd, order, ii);
    end

    tbx_sm_plot_simu(res, is_saved);
