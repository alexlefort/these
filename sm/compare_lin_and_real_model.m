
Vs = 1.5;

param   = tbx_sm_load_model('iver2.nav') ;
noise   = tbx_sm_load_noises()           ;
mission = tbx_sm_load_mission(Vs)        ;
param_lin = param;
dt = 0.1;

    rng(1234);
    close all;
    
    %% Initial State
    
    state  = mission.state_init ;
    cmd    = mission.cmd_init   ;
    res    = tbx_sm_make_simu_vect(mission.T);
    state_lin = state;
    
    %% Simulation Loop
    t = 0.0;
    
    for ii=1:(length(mission.T))
         t = t + dt;
         
         cmd.Alpha = 0;
         if (t > 1) 
             cmd.Alpha = 0.1;
         end
         state.U = Vs;
         state        = tbx_sm_dynamics(state, cmd, param, dt);
         state_lin    = tbx_sm_dynamics_lin_gir(state_lin, cmd, param_lin, dt);
         measure      = tbx_sm_sensors(state, noise,0);        
         res          = tbx_sm_save_simu_vect(res, state, measure, cmd, order, ii);
         res_lin      = tbx_sm_save_simu_vect(res_lin, state_lin, measure, cmd, order, ii);
    end
    
    T = mission.T;
    R = res.state.R;
    R_lin = res_lin.state.R;
    figure;
    plot(T,R,'r',T,R_lin,'b');