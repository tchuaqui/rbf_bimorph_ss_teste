clear all

L=100e-3;
ni=21;    %DEFINICAO DO NUMERO DE NOS    (n=200 com c=2/sqrt(n) É O QUE DÁ MELHOR!!!!)
x=0:L/(ni-1):L;n=numel(x);

cfstr='ss'; %DEFINICAO DAS CONDICOES DE FRONTEIRA - CC OU SS
d=1;                     %DEFINICAO DO DIAMETRO OU DA ESPESSURA DA VIGA
carga=1000;                          %DEFINICAO DA CARGA
% A=b*h;               
k=5/6;
%propriedades do piezoelectrico
propmec=[5e-3 6e10 2.3e10 7500]; %[esp E G rho]
propel=[-16.492145 -16.492145 2.588549e-8]; %[e31 e32 ezz]
esp=[propmec(1) propmec(1)];
rho=propmec(4);

%%
e31=propel(1);
ezz=propel(3);

Q11=[propmec(2) propmec(2)];
Q55=[propmec(3) propmec(3)];
%%
ncamadas=2;
h=0;
for i=1:ncamadas
h=h+esp(i);    
end

z=zeros(ncamadas+1,1);
z(1)=-h/2;
zm=zeros(ncamadas,1);
for i=2:ncamadas+1
z(i)=z(i-1)+esp(i-1);
zm(i-1)=(z(i-1)+z(i))/2;
end

I0=zeros(ncamadas,1); I1=zeros(ncamadas,1); I2=zeros(ncamadas,1); J0=0; J1=0; J2=0;
B11=0; B55=0; C11=0; D11=0; 
for i=1:ncamadas    
I0(i)=z(i+1)-z(i);
I1(i)=(z(i+1)^2-z(i)^2)/2;
I2(i)=(z(i+1)^3-z(i)^3)/3;

J0=J0+rho*I0(i);
J1=J1+rho*I1(i);
J2=J2+rho*I2(i);

B11=B11+Q11(i)*I0(i);
C11=C11+Q11(i)*I1(i);
B55=B55+k*Q55(i)*I0(i);
% Atheta=Atheta+Q11(i)*I1(i);
% Btheta=Btheta-k*Q55(i)*I0(i);
D11=D11+Q11(i)*I2(i);   %c1theta
% C2theta=C2theta-k*Q55(i)*I0(i); 
% Mu=Mu+Q11(i)*I1(i);
% Mtheta=Mtheta+Q11(i)*I2(i);
end

F1=e31*I0(1)/esp(1);
F2=e31*I0(2)/esp(2);
H1=e31*I1(1)/esp(1);
H2=e31*I1(2)/esp(2);
I_1=ezz*I0(1)/(esp(1)^2);
I_2=ezz*I0(2)/(esp(2)^2);

%%
c=L./sqrt(numel(x)); 
% c=2./sqrt(numel(x)); 
% c=2./numel(x);
% c=cross_validation(cfstr,n,L,x,Au,Cu,Bw,Atheta,Btheta,C1theta,C2theta,vetor_carga,Du1,Du2,Dtheta1,Dtheta2,phi1,phi2,Mu,Mtheta);
%%

[xi,xj]=meshgrid(x,x);
AA=g(c,xi,xj);
DAA=dgdx(c,xi,xj);
D2AA=d2gdx2(c,xi,xj);


digits(500)

S_total(1:n,1:n)=B11*D2AA;
S_total(1:n,n+1:2*n)=0;
S_total(1:n,2*n+1:3*n)=C11*D2AA;

A_total(1:n,1:n)=-J0*AA;
A_total(1:n,n+1:2*n)=0;
A_total(1:n,2*n+1:3*n)=-J1*AA;

S_total(n+1:2*n,1:n)=0;                    
S_total(n+1:2*n,n+1:2*n)=B55*D2AA;
S_total(n+1:2*n,2*n+1:3*n)=B55*DAA;

A_total(n+1:2*n,1:n)=0;                    
A_total(n+1:2*n,n+1:2*n)=-J0*AA;
A_total(n+1:2*n,2*n+1:3*n)=0;

S_total(2*n+1:3*n,1:n)=C11*D2AA;
S_total(2*n+1:3*n,n+1:2*n)=-B55*DAA;
S_total(2*n+1:3*n,2*n+1:3*n)=D11*D2AA-B55*AA;

A_total(2*n+1:3*n,1:n)=-J1*AA;
A_total(2*n+1:3*n,n+1:2*n)=0;
A_total(2*n+1:3*n,2*n+1:3*n)=-J2*AA;

Kuphis(1:n,1:n)=F1*DAA;
Kuphia(1:n,1:n)=F2*DAA;

Ktphis(1:n,1:n)=H1*DAA;
Ktphia(1:n,1:n)=H2*DAA;   

Kphiphis(1:n,1:n)=I_1*AA;
Kphiphia(1:n,1:n)=I_2*AA;


Kuu=S_total(1:n,1:n)+Kuphis*(Kphiphis^-1)*Kuphis;
Kut=S_total(1:n,2*n+1:3*n)+Kuphis*(Kphiphis^-1)*Ktphis;
Kww=S_total(n+1:2*n,n+1:2*n);
Kwt=S_total(n+1:2*n,2*n+1:3*n);
Ktu=S_total(2*n+1:3*n,1:n)+Ktphis*(Kphiphis^-1)*Kuphis;
Ktw=S_total(2*n+1:3*n,n+1:2*n);
Ktt=S_total(2*n+1:3*n,2*n+1:3*n)+Ktphis*(Kphiphis^-1)*Ktphis;

K_total(1:n,1:n)=Kuu; K_total(1:n,n+1:2*n)=0; K_total(1:n,2*n+1:3*n)=Kut;
K_total(n+1:2*n,1:n)=0; K_total(n+1:2*n,n+1:2*n)=Kww; K_total(n+1:2*n,2*n+1:3*n)=Kwt;
K_total(2*n+1:3*n,1:n)=Ktu; K_total(2*n+1:3*n,n+1:2*n)=Ktw; K_total(2*n+1:3*n,2*n+1:3*n)=Ktt;


%CONDICOES DE FRONTEIRA

switch cfstr
    
    case {'ss'}
    
%equacao de fronteira - CASO SIMPLESMENTE APOIADO
b=find(x==0 | x==L );
% S_total(b,:)=[Au*DAA(b,:),zeros(2,n),Cu*DAA(b,:)]; %PARA NX
K_total(b,:)=[AA(b,:),zeros(2,n),zeros(2,n)];         %PARA U
K_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
K_total(b+2*n,:)=[C11*DAA(b,:),zeros(2,n),D11*DAA(b,:)]; %PARA M

A_total(b,:)=zeros(2,3*n);
A_total(b+n,:)=zeros(2,3*n);
A_total(b+2*n,:)=zeros(2,3*n);

    case {'cc'}        
%equacao de fronteira - CASO ENCASTRADO
b=find(x==0 | x==L );
K_total(b,:)=[AA(b,:),zeros(2,n),zeros(2,n)];          %PARA U
K_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
K_total(b+2*n,:)=[zeros(2,n),zeros(2,n),AA(b,:)];                 %PARA THETA

A_total(b,:)=zeros(2,3*n);
A_total(b+n,:)=zeros(2,3*n);
A_total(b+2*n,:)=zeros(2,3*n);

    case {'cl'}
b1=find(x==0); 
K_total(b1,:)=[AA(b1,:),zeros(1,n),zeros(1,n)];          %PARA U
K_total(b1+n,:)=[zeros(1,n),AA(b1,:),zeros(1,n)];       %PARA W
K_total(b1+2*n,:)=[zeros(1,n),zeros(1,n),AA(b1,:)];                 %PARA THETA

A_total(b1,:)=zeros(1,3*n);
A_total(b1+n,:)=zeros(1,3*n);
A_total(b1+2*n,:)=zeros(1,3*n);

b2=find(x==L);
K_total(b2,:)=[Au*DAA(b2,:),zeros(1,n),Cu*DAA(b2,:)]; %PARA NX
K_total(b2+n,:)=[zeros(1,n),Bw*DAA(b2,:),Bw*AA(b2,:)];       %PARA QX
K_total(b2+2*n,:)=[Mu*DAA(b2,:),zeros(1,n),Mtheta*DAA(b2,:)];       %PARA MX

A_total(b2,:)=zeros(1,3*n);
A_total(b2+n,:)=zeros(1,3*n);
A_total(b2+2*n,:)=zeros(1,3*n);
end


%% DINAMICA

[lambda_vec,lambda]=eig(K_total,A_total);

lambda=diag(lambda,0); 
[lambda,indice]=sort(lambda);
lambda_vec=lambda_vec(:,(indice(:)));
lambda_norm=sqrt(lambda)*L^2*sqrt(rho*h/(Q11*I2));

% lambda=real(lambda);
%  lambda(lambda(:)==-inf | lambda(:)<=0)=NaN;
% 
% lambda=sqrt(lambda)*L^2*sqrt(rho*A/(E*I));%lambda=sqrt(lambda);


m=1; E=6e10; I=I2(1)+I2(2); A=h; G=2.3e10;
sol_exacta=(m*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((m*pi/L)^2*E*I)/(k*G*A+(m*pi/L)^2*E*I)));
% sol_exacta_norm=sol_exacta*L^2*sqrt(rho*A/(E*I));
%-------------graficos 


p=1;
lambda_mode_w(1:n,p)=lambda_vec(n+1:2*n,p)'*AA;
lambda_mode_phi_x(1:n,p)=lambda_vec(2*n+1:3*n,p)'*AA;
lambda_mode=[lambda_mode_w;lambda_mode_phi_x];


% SS3=lambda_vec(:,1)'*A_total*lambda_vec(:,1);
% % SS=*lambda_vec(:,p)-lambda(p,p)*A_total*lambda_vec(:,p);
% SS=(S_total+lambda(1)*A_total)*lambda_vec(:,1);
% 
% SS4=S_total*real(lambda_vec3)-A_total*real(lambda_vec3)*real(lambda3);
% SS4=S_total*(lambda_vec3)-A_total*(lambda3);

% subplot(1,3,1);plot(x, lambda_mode_w(:,p));hold on;
% % title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta_norm,'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(lambda(m),'%6.4f')])
% subplot(1,3,2);plot(x, lambda_mode_phi_x(:,p));hold on;
% % title(['w(' num2str(m) ')_{exact} = ' num2str(sol_exacta_norm,'%6.4f')]);legend(['w(' num2str(m) ')=' num2str(lambda(m),'%6.4f')])
freq=sqrt(lambda(p))/(2*pi);

figure(1)
plot(x,lambda_mode_w);
%NOTA: ERRO CALCULADO ABAIXO CORRESPONDE AO ERRO RELATIVO COM OS RESULTADOS
%DO REDDY PARA UMA VIGA SS NAO PIEZO
plot(x, lambda_mode_w(:,p));hold on;title(['Forma natural' num2str(p)]);legend(['wnorm(' num2str(p) ')=' num2str(lambda_norm(p),'%6.4f') ', w(' num2str(p) ')=' num2str(sqrt(lambda(p))/(2*pi),'%6.4f') 'Hz' ', w exata=' num2str(sol_exacta/(2*pi),'%6.4f') 'Hz' ,', erro=' num2str(100*(sol_exacta-sqrt(lambda(p)))/(sol_exacta),'%6.4f') '%'])
%%
%-----NEWMARK

%CONTROLADOR
Gv=0.001;

Cuu=Kuphia*(Kphiphis^-1)*Kuphis;
Cut=Kuphia*(Kphiphis^-1)*Ktphis;
Ctu=Ktphia*(Kphiphis^-1)*Kuphis;
Ctt=Ktphia*(Kphiphis^-1)*Ktphis;

C_total=zeros(3*n,3*n);
C_total(1:n,1:n)=-Cuu;   %estava positivo, com negativo funciona!!!!!!!!!!!!!!!!!
C_total(1:n,2*n+1:3*n)=Cut;
C_total(2*n+1:3*n,1:n)=Ctu;
C_total(2*n+1:3*n,2*n+1:3*n)=Ctt;
C_total=-Gv*C_total;

b=find(x==0 | x==L );
C_total(b,:)=zeros(2,3*n);
C_total(b+n,:)=zeros(2,3*n);
C_total(b+2*n,:)=zeros(2,3*n);
%

%cond. iniciais 
vetor_carga=zeros(3*n,1);
vetor_carga(n+2:2*n-1)=carga;
solucao_estatica=K_total\vetor_carga;
x_0=solucao_estatica; v_0=zeros(3*n,1);
vetor_f=zeros(3*n,1);
a_0=pinv(A_total)*(vetor_f-C_total*v_0-K_total*x_0);
delta=1/2; alpha=1/4;
% delta_t=1/(freq*200);   %delta t
delta_t=1/100000;
a0=1/(alpha*delta_t^2); a1=delta/(alpha*delta_t); a2=1/(alpha*delta_t); a3=1/(2*alpha)-1;
a4=delta/alpha-1; a5=(delta_t/2)*(delta/alpha-2); a6=delta_t*(1-delta); a7=delta*delta_t;

%rigidez efetiva
K_efe=K_total+a0*A_total+a1*C_total;

% t_final=10/freq;   % t final
t_final=0.02;
n_t=int64(t_final/delta_t+1);
t=zeros(n_t,1);
x_t=zeros(3*n,n_t); x_t(:,1)=x_0;
v_t=zeros(3*n,n_t); v_t(:,1)=v_0;
a_t=zeros(3*n,n_t); a_t(:,1)=a_0;
for i=2:n_t
  t(i)=t(i-1)+delta_t;
  F_efe=A_total*(a0*x_t(:,i-1)+a2*v_t(:,i-1)+a3*a_t(:,i-1))+C_total*(a1*x_t(:,i-1)+a4*v_t(:,i-1)+a5*a_t(:,i-1));
  x_t(:,i)=K_efe\F_efe;
  a_t(:,i)=a0*(x_t(:,i)-x_t(:,i-1))-a2*v_t(:,i-1)-a3*a_t(:,i-1);
  v_t(:,i)=v_t(:,i-1)+a6*a_t(:,i-1)+a7*a_t(:,i);
end

%interpolacao
x_w=zeros(n,n_t);
x_max=zeros(n_t,1);
for i=1:n_t
aux=x_t(n+1:2*n,i)'*AA;
x_w(:,i)=aux';
x_max(i)=x_w(ceil(n/2),i);
end

figure(2)
plot(t,x_max);
hold on

X=fft(x_max);
X_mag=abs(X(1:ceil(n_t/2)));
[pk_vals, pk_locs]=findpeaks(X_mag);
%remove peaks below threshold
% inds=find(X_mag(pk_locs)<1);
% pk_locs(inds)=[];

%determine frequencies
pk_freqs=zeros(length(pk_locs),1);
for i=1:length(pk_locs)
pk_freqs(i)=(pk_locs(i)-1)/t_final;
end

figure(3)
plot(X_mag);
hold on
