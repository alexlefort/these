%% Return the hinf criteria of a problem as its modulus
%%

function criteria_hinf = build_criteria_hinf(mod,p,k)

    %% Define all parameters as real values

	params = fields(p);
	nparams = length(fields(p));

	for kkk=1:nparams
		scrip = strcat(params{kkk}, ' = sym(''',params{kkk},''',''real'');');
		eval(scrip);
		mod = subs(mod, params{kkk}, params{kkk});
	end

    %% Define all gains as real values

	gains = fields(k);
	ngains = length(fields(k));

	for kkk=1:ngains
		scrip = strcat(gains{kkk}, ' = sym(''',gains{kkk},''',''real'');');
		eval(scrip);
		mod = subs(mod, gains{kkk}, gains{kkk});
	end

    %% Change s laplace variable to pulsation 1i*omega (w)

	w = sym('w','real');
	mod = subs(mod, 's', 1i*w);
	
	%% Square module of the transfert function as hinf criteria

	[num_mod,den_mod] = numden(mod);

    numD = real(num_mod)^2 + imag(num_mod)^2;
    denD = real(den_mod)^2 + imag(den_mod)^2;

    criteria_hinf = numD/denD;