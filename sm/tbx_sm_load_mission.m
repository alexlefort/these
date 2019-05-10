function mission = tbx_sm_load_mission()

    dt        = 0.1        ;
    Tmax      = 1000       ;
    mission.T = 0:dt:Tmax  ;
    
    n = length(mission.T);
    
    mission.state_init.U     = 0.5 ;
    mission.state_init.V     = 0.0 ;
    mission.state_init.W     = 0.0 ;
    mission.state_init.P     = 0.0 ;
    mission.state_init.Q     = 0.0 ;
    mission.state_init.R     = 0.0 ;
    mission.state_init.X     = 0.0 ;
    mission.state_init.Y     = 0.0 ;
    mission.state_init.Z     = 0.0 ;
    mission.state_init.Phi   = 0.0 ; 
    mission.state_init.Theta = 0.0 ;
    mission.state_init.Psi   = 0.0 ; 
    
    mission.cmd_init.Alpha  =  0.0 ;
    mission.cmd_init.Beta1  =  0.0 ;
    mission.cmd_init.Beta2  =  0.0 ;
    mission.cmd_init.Prop   = 12.0 ;
   
    mission.ZCo   = 50*ones(n,1) ;
    mission.PsiCo = 1*ones(n,1) ;
    mission.Prop  = 12*ones(n,1);