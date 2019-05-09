function measure = tbx_sm_sensors(state, noise)

    measure.U     = state.U*randn(noise.U);
    measure.V     = state.V*randn(noise.V);
    measure.P     = state.P*randn(noise.P);
    measure.Q     = state.Q*randn(noise.Q);
    measure.R     = state.R*randn(noise.R);
    measure.Phi   = state.Phi*randn(noise.Phi);
    measure.Theta = state.Theta*randn(noise.Theta);
    measure.Psi   = state.Psi*randn(noise.Psi);
    measure.Z     = state.Z*randn(noise.Z);