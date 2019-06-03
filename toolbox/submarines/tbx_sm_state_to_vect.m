function xvect = tbx_sm_state_to_vect(xstate)

    xvect(1,1)  = xstate.U     ;
    xvect(2,1)  = xstate.V     ;
    xvect(3,1)  = xstate.W     ;
    xvect(4,1)  = xstate.P     ;
    xvect(5,1)  = xstate.Q     ;
    xvect(6,1)  = xstate.R     ;
    xvect(7,1)  = xstate.X     ;
    xvect(8,1)  = xstate.Y     ;
    xvect(9,1)  = xstate.Z     ;
    xvect(10,1) = xstate.Phi   ;
    xvect(11,1) = xstate.Theta ;
    xvect(12,1) = xstate.Psi   ;