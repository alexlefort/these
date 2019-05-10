function [X,U] = init_etats
states = {'U';'V';'W';'P';'Q';'R';'Phi';'Theta';'Psi';'Z';'fxp';'fyp';'fzp';'mxp';'myp';'mzp'};
inputs = {'Alpha';'Beta1';'Beta2';'NT1'};

for jj=1:numel(states)
    X.(states{jj}) = 0.0;
end

for jj=1:numel(inputs)
    U.(inputs{jj}) = 0.0;
end
