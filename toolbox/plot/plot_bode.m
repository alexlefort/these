%% Plot a bode diagram for a set of uncertain parameter for a
%% Symbolic transfer function

function plot_bode(sym_transfer, param_vect, logw, sizefont, legend_t, namefile)

    syms s;

    w = 10.^(logw); %% Pulsation vector
    
    param_names = fields(param_vect);
    param_table = param_org_table(param_vect);
    [nb_points, nb_param] = size(param_table);
    
    param_val = param_vect;
    
    for ii=1:nb_param
        param_val.(param_names{ii}) = 0;
    end

    fg = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    x   = xlabel('Pulsation (rad/s)');
    set(x,'Interpreter','latex','FontSize',sizefont);
    y   = ylabel('Magnitude (dB)');
    set(y,'Interpreter','latex','FontSize',sizefont);

    %% Plot closed-loop transfer functions

    for ii=1:nb_points
        for jj=1:nb_param
            param_val.(param_names{jj}) = param_table(ii,jj);
        end
        
        transfer_0 = symtbx_sym_subs_from_struct(sym_transfer, param_val); 
        bode_0     = symtbx_symtf_bode(transfer_0, w); 

        hold on;
        set(gca, 'YScale', 'log')
        set(gca, 'XScale', 'log')
        loglog(w, bode_0, 'b'); grid on;

    end
    h = legend(legend_t);
    set(h,'Interpreter','latex','FontSize',sizefont);
    
    if (nargin > 6)
        saveas(fg, namefile);
    end