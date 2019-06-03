function simulation(sys,ctrl,T)
    
    deg = pi/180;
    assignin('base', 'sys'  , sys  );
    assignin('base', 'ctrl' , ctrl );

    sim('ctrl_imm_simu', 0:ctrl.dt:(T));

    sorties = { 'Z' ;'Theta' ; 'Pi' ; 'Q' ; 'Beta'};
    unites  = {  1  ; deg    ;  deg ; deg ;   deg };
    nb_sorties = length(sorties);
    tt = output.Z.Time;
    
    figure;
    for kk=1:nb_sorties
        sortie = sorties{kk};
        subplot(nb_sorties,1,kk);
        data = output.(sortie).Data/unites{kk};
        plot(tt,data,'b','linewidth',3);
        ylabel(sortie);
        grid on;
    end
    hold off;
