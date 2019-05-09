function res = simulation()
    

    param      = get_model()      ;
    noise      = get_noise()      ;
    mission    = get_mission()    ;
    param_ctrl = get_param_ctrl() ;
    
    %% Initial State
    state_init   = mission.state_init;
    command_init = mission.command_init;

    %% Simulation Loop
    
    state   = state_init   ;
    command = command_init ;
    
    res = make_vect(T);

    for ii=1:length(T)
         t = t + param.dt;
         state   = sm_dynamics(state, command, param, param.dt);
         measure = sensors(state, noise, param.dt);
         command = controller(measure, param_ctrl, intern, param.dt);
         res = save_vect(res, state, measure, command, mission, T);
    end


