function res = tbx_sm_simulation(is_saved)
    
    rng(1234);
   
    param   = tbx_sm_load_model()   ;
    noise   = tbx_sm_load_noise()   ;
    mission = tbx_sm_load_mission() ;
    pctrl   = tbx_sm_load_pctrl()   ;
    
    %% Initial State
    
    state   = mission.state_init   ;
    command = mission.command_init ;
    intern  = pctrl.intern_init    ;
    
    res = tbx_sm_make_simu_vect(T);
    
    %% Simulation Loop

    for ii=1:length(T)
         t = t + param.dt;
         state   = tbx_sm_dynamics(state, command, param, param.dt);
         measure = tbx_sm_sensors(state, noise);
         command = tbx_sm_controller(measure, pctrl, intern, param.dt);
         res     = tbx_sm_save_simu_vect(res, state, measure, command, mission, ii);
    end

    tbx_sm_plot_simu(res, is_saved);
