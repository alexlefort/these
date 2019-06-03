function measure = tbx_sm_sensors(state, noise, active)

    measure.U     = state.U*(1+randn(1)*noise.U*active);
    measure.V     = state.V*(1+randn(1)*noise.V*active);
    measure.W     = state.W*(1+randn(1)*noise.W*active);
    measure.P     = state.P*(1+randn(1)*noise.P*active);
    measure.Q     = state.Q*(1+randn(1)*noise.Q*active);
    measure.R     = state.R*(1+randn(1)*noise.R*active);
    measure.Phi   = state.Phi*(1+randn(1)*noise.Phi*active);
    measure.Theta = state.Theta*(1+randn(1)*noise.Theta*active);
    measure.Psi   = state.Psi*(1+randn(1)*noise.Psi*active);
    measure.Z     = state.Z*(1+randn(1)*noise.Z*active);