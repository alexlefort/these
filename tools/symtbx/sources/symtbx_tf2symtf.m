%% Transform an expression into a matlab transfer function
%%

function symtf = symtbx_tf2symtf(tf_t,s)

   [num,den] = tfdata(tf_t);
   symtf = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);