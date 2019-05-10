function measure = tbx_sm_sensors(state, noise)

    measure.U     = state.U*(1+randn(1)*noise.U);
    measure.V     = state.V*(1+randn(1)*noise.V);
    measure.W     = state.W*(1+randn(1)*noise.W);
    measure.P     = state.P*(1+randn(1)*noise.P);
    measure.Q     = state.Q*(1+randn(1)*noise.Q);
    measure.R     = state.R*(1+randn(1)*noise.R);
    measure.Phi   = state.Phi*(1+randn(1)*noise.Phi);
    measure.Theta = state.Theta*(1+randn(1)*noise.Theta);
    measure.Psi   = state.Psi*(1+randn(1)*noise.Psi);
    measure.Z     = state.Z*(1+randn(1)*noise.Z);