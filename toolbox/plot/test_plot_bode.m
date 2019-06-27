%% Test plot_bode function

clear all
close all

syms s p1 p2 p3;

p.p1 = -0.235:0.03:0.235;
p.p2 = 0.998:0.001:1.002;
p.p3 = 7;

sym_transfer = (p3*(s-2))/(s^2 + p2*s + p1);
logw = -3:0.01:1;
sizefont = 14;
legend_t = '$T_{b \rightarrow z}(p^*)$';
namefile = 'test_plot_bode';

plot_bode(sym_transfer, p, logw, sizefont, legend_t, namefile);
