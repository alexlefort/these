f = model_capd();

p.CmQ    = -0.8116   ;
p.CzW    = -2.6126   ;
p.CmB1   = -0.7155   ;
p.CzB1   = -1.2233   ;
p.ZCo    =  0.0      ;
p.kpi    = -0.033225 ;
p.kq     = -0.33076  ;
p.ktheta = -0.82888  ;
p.kz     = -0.1709  ;
p.iz     = 0.0 ;

x.Z     = 0;
x.Theta = 0;
x.Pi    = 0;
x.Q     = 0;
x.IZ    = 0;
x.Beta  = 0;

f0 =  symtbx_sym_subs_from_struct(f, p);

n = length(f0);

A = zeros(n);

for ii=1:length(f0)
    z = zeros(6,1);
    z(ii) = 1;
    x.Z     = z(1);
    x.Theta = z(2);
    x.Pi    = z(3);
    x.Q     = z(4);
    x.IZ    = z(5);
    x.Beta  = z(6);
    
    for jj=1:length(f0)
        res = eval(symtbx_sym_subs_from_struct(f0, x));
        A(jj,ii) = res(jj);
    end
end

tspan = 1:0.01:10;
x0 = [1;0;0;0;0;0];

[t,y] = ode45(@myode,tspan,x0);
figure;
subplot(411);plot(t,y(:,1));
subplot(412);plot(t,y(:,2));
subplot(413);plot(t,y(:,3));
subplot(414);plot(t,y(:,4));

