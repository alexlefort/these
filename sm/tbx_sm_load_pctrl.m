function pctrl = tbx_sm_load_pctrl()

    aux1 = load('param_depth_controller.mat'  );
    aux2 = load('param_heading_controller.mat');
    
    pctrl.depth = aux1.param_depth_controller;
    pctrl.heading = aux2.param_heading_controller;
    
    pctrl.dt                         = 0.1 ;
    pctrl.intern_init.heading.Alpha  = 0.0 ;
    pctrl.intern_init.heading.Rint   = 0.0 ;
    pctrl.intern_init.heading.antiwp = 0.0 ;
    
    pctrl.intern_init.depth.Beta   = 0.0 ;
    pctrl.intern_init.depth.Zint   = 0.0 ;
    pctrl.intern_init.depth.antiwp = 0.0 ;