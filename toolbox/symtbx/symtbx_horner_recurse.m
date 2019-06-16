%% Recurse function for horner algorithm
%%

function p_new = symtbx_horner_recurse(p)

    val   = symvar(p);
    nval  = length(val);

    if nval>=1
        for ii = 1:nval
            [monomials{ii}, orders{ii}] = coeffs(p,val(ii));
            noccurs(ii) = length(orders{ii});
            strr = char(monomials{ii}(1));
            nsums(ii) = length(strfind(strr,'+'));
        end

        [valoccur, idx] = min(noccurs);
        pivot = val(idx);

        if (valoccur > 1)
            [aux,idx2] = max(nsums);
            pivot = val(idx2);
        end

        if (valoccur >= 1)
    
            [new_monomial, new_order] = coeffs(p,pivot);

            %% Decomposition p = p1*pivot + p2;
    
            p1 = 0;
            for jj=1:(valoccur)
    
                if (new_order(jj) ~= 1)
                    p1 = p1 + new_monomial(jj)*new_order(jj)/pivot;
                end
            end
    
            p1_new = symtbx_horner_recurse(p1);
            p_new = p1_new*pivot;
    
            if (new_order(valoccur) == 1)
                p2 = new_monomial(valoccur);
                p2 = symtbx_horner_recurse(p2);
                p_new = p_new + p2;
            end
    
        else
            p_new = p;
        end
    else
        p_new = simplify(p);
    end
