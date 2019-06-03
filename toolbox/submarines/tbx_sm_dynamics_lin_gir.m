function x = tbx_sm_dynamics_lin_gir(x,u,p,dt)

    x0 = tbx_sm_state_to_vect(x);
    k1 = model_lin_gir(x0, u, p);    
    k2 = model_lin_gir(x0 + dt/2*k1 , u, p);
    k3 = model_lin_gir(x0 + dt/2*k2 , u, p);
    k4 = model_lin_gir(x0 + dt*k3   , u, p);
    x  = x0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    x  = tbx_sm_vect_to_state(x);


function xpt = model_lin_gir(x,u,p)
    xpt = zeros(12,1);
    p.Vs = x(1);
    x_gir = x([12 2 6]);
    u_gir = u.Alpha;
    model = tbx_sm_heading_model(p);
    xpt_gir = model.a*x_gir + model.b*u_gir;
    xpt(12) = xpt_gir(1);
    xpt(2)  = xpt_gir(2);
    xpt(6)  = xpt_gir(3);