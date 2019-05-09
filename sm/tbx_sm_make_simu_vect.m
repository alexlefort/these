function res = tbx_sm_make_simu_vect(T)

    n = length(T);
    
    res.state.U       = zeros(n,1);
    res.state.V       = zeros(n,1);
    res.state.W       = zeros(n,1);
    res.state.P       = zeros(n,1);
    res.state.Q       = zeros(n,1);
    res.state.R       = zeros(n,1);
    res.state.X       = zeros(n,1);
    res.state.Y       = zeros(n,1);
    res.state.Z       = zeros(n,1);
    res.state.Phi     = zeros(n,1);
    res.state.Theta   = zeros(n,1);
    res.state.Psi     = zeros(n,1);
    
    res.measure.U     = zeros(n,1);
    res.measure.V     = zeros(n,1);
    res.measure.P     = zeros(n,1);
    res.measure.Q     = zeros(n,1);
    res.measure.R     = zeros(n,1);
    res.measure.Z     = zeros(n,1);
    res.measure.Phi   = zeros(n,1);
    res.measure.Theta = zeros(n,1);
    res.measure.Psi   = zeros(n,1);
    
    res.command.Alpha  = zeros(n,1);
    res.command.Beta1  = zeros(n,1);
    res.command.Beta2  = zeros(n,1);
    res.command.Prop   = zeros(n,1);
    
    res.mission.ZCo    = zeros(n,1);
    res.mission.PsiCo  = zeros(n,1);
    res.mission.PropCo = zeros(n,1);
    