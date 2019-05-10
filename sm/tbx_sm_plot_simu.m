function tbx_sm_plot_simu(res, is_saved)

figure; plot(res.T,res.state.Z);
figure; plot(res.T,res.state.Psi);