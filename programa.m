clear all
L=100e-3; nmax=100;
omega=zeros(nmax-1,1); omega_o=zeros(nmax-1,1);
p=1;

for ni=2:nmax
 x=0:L/(ni-1):L; n=numel(x);
%  
c=L./sqrt(numel(x)); 
% c=2./sqrt(numel(x)); 
% c=2./(numel(x));
% c=2/L;

 [omega(ni-1),E,I,A,G,rho,k] = timo_linear_clo( ni,L,c,p);   
[omega_o(ni-1),E,I,A,G,rho,k] = timo_linear_op( ni,L,c,p); 
end

sol_exacta=(p*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((p*pi/L)^2*E*I)/(k*G*A+(p*pi/L)^2*E*I)));
sol_exacta=sol_exacta/(2*pi); 
sol_ex_vec=sol_exacta*ones(length(omega)+2,1);
figure(1)
plot(2:nmax,omega);
hold on
plot(2:nmax,omega_o);
hold on
plot(0:nmax,sol_ex_vec);
