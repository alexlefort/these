function xvect = tbx_sm_state_to_vect(xstate)

    xvect(1)  = xstate.U     ;
    xvect(2)  = xstate.V     ;
    xvect(3)  = xstate.W     ;
    xvect(4)  = xstate.P     ;
    xvect(5)  = xstate.Q     ;
    xvect(6)  = xstate.R     ;
    xvect(7)  = xstate.X     ;
    xvect(8)  = xstate.Y     ;
    xvect(9)  = xstate.Z     ;
    xvect(10) = xstate.Phi   ;
    xvect(11) = xstate.Theta ;
    xvect(12) = xstate.Psi   ;