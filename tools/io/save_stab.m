function save_stab(c, z, k, p,name)


	file = fopen(name,'w');
	
	kv  = fields(k);
	nkv = size(kv,1);

	pv  = fields(p);
	npv = size(pv,1);

	for ii=1:nkv
		ii
		kv{ii}
	end 

for ii=1:npv
		ii
		pv{ii}
	end 

	%% Write head of function
	
	mes = 'function f(';
		
	mes = strcat(mes , 'x[' , num2str(nkv) , '], p[' , num2str(npv) , ']) \n return(');

	%% Write all polys

	pv  = fields(p);
	npv = size(pv,1);
		
	n = size(c,2);
	
	for ii = 1:n
		expr = char(c{ii});
			
		%% Substitute generic names for parameters and gains (xi,pi) :
		
		for jj=1:nkv
			char_x = strcat('x[',num2str(jj-1),']');
			expr   = strrep(expr ,kv{jj}, char_x);
		end

		for jj=1:npv
			char_p = strcat('p[',num2str(jj-1),']');
			expr   = strrep(expr ,pv{jj}, char_p);
		end
		
		mes = strcat(mes,expr);
        if ii<n
            mes = strcat(mes,',');
        end
	
	end

    mes = strcat(mes,') \n end');
			
	fprintf(file,mes);
	fclose(file);
