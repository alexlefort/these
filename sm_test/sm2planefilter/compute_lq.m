clear all
clc
close all


model_sm = build_model_lin();

ctrl_gains.kpt1 = -0.899498;
ctrl_gains.kdt1 = -6.94992;
ctrl_gains.kiz2 =  0.019999;
ctrl_gains.kpz2 =  0.475193;
ctrl_gains.kdz2 = -2.34777;


syms s;

logw = -2:0.1:1;
w    = 10.^(logw); %% Pulsation vector

 disp(ctrl_gains);
 
 load('criterias.mat');
 
 czz1 = symtbx_sym_subs_from_struct(criterias.zz1, model_sm.p0);
 czz1 = simplify(symtbx_sym_subs_from_struct(czz1, ctrl_gains));

 czz1_bode  = eval(max(abs(subs(czz1,'w', w))));

 disp(czz1_bode);
 
 czz2 = symtbx_sym_subs_from_struct(criterias.zz2, model_sm.p0);
 czz2 = simplify(symtbx_sym_subs_from_struct(czz2, ctrl_gains));
 czz2_bode  = eval(max(abs(subs(czz2,'w', w))));
 
 disp(czz2_bode);
 
 czb1 = symtbx_sym_subs_from_struct(criterias.zb1, model_sm.p0);
 czb1 = simplify(symtbx_sym_subs_from_struct(czb1, ctrl_gains));
 czb1_bode  = eval(max(abs(subs(czb1, 'w', w))));
 
 disp(czb1_bode);
 
 czb2 = symtbx_sym_subs_from_struct(criterias.zb2, model_sm.p0);
 czb2 = simplify(symtbx_sym_subs_from_struct(czb2, ctrl_gains));
 czb2_bode  = eval(max(abs(subs(czb2, 'w', w))));
 
 disp(czb1_bode);

for ii=1:5
        crit = criterias.stab_lc{ii};
        crit = symtbx_sym_subs_from_struct(crit, model_sm.p0);
        crit = symtbx_sym_subs_from_struct(crit, ctrl_gains);
        disp(eval(crit));
end

disp(max(max(max(max(czz1_bode),max(czz2_bode)),max(max(czb1_bode),max(czb1_bode)))));
