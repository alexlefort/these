function [Asub,Bsub] = symauto_blockmod(A,B,all,sub)


	Asub = sym(zeros(sub.nb_states, sub.nb_states));
	
	for ii=1:sub.nb_states
		for jj=1:sub.nb_states
			for kk=1:all.nb_states
				for ll=1:all.nb_states
					if (strcmp(sub.states{ii},all.states{kk}) && strcmp(sub.states{jj},all.states{ll}))
						Asub(ii,jj) = A(kk,ll);
					end
				end
			end
		end
	end
	
	Bsub = sym(zeros(sub.nb_states, sub.nb_inputs));
	
	for ii=1:sub.nb_states
		for jj=1:sub.nb_inputs
			for kk=1:all.nb_states
				for ll=1:all.nb_inputs
					if (strcmp(sub.states{ii},all.states{kk}) && strcmp(sub.inputs{jj},all.inputs{ll}))
						Bsub(ii,jj) = B(kk,ll);
					end
				end
			end
		end
	end
