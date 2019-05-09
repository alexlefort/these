function xstate = tbx_sm_vect_to_state(xvect)

    xstate.U      = xvect(1)  ;
    xstate.V      = xvect(2)  ;
    xstate.W      = xvect(3)  ;
    xstate.P      = xvect(4)  ;
    xstate.Q      = xvect(5)  ;
    xstate.R      = xvect(6)  ;
    xstate.X      = xvect(7)  ;
    xstate.Y      = xvect(8)  ;
    xstate.Z      = xvect(9)  ;
    xstate.Phi    = xvect(10) ;
    xstate.Theta  = xvect(11) ;
    xstate.Psi    = xvect(12) ;