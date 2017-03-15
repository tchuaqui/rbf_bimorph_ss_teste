clear all

L=100e-3;
ni=30;    %DEFINICAO DO NUMERO DE NOS    (n=200 com c=2/sqrt(n) É O QUE DÁ MELHOR!!!!)
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
Au=0; Cu=0; Bw=0; Atheta=0; Btheta=0; C1theta=0; C2theta=0; Mu=0; Mtheta=0;
for i=1:ncamadas    
I0(i)=z(i+1)-z(i);
I1(i)=(z(i+1)^2-z(i)^2)/2;
I2(i)=(z(i+1)^3-z(i)^3)/3;

J0=J0-rho*I0(i);
J1=J1-rho*I1(i);
J2=J2-rho*I2(i);

Au=Au+Q11(i)*I0(i);
Cu=Cu+Q11(i)*I1(i);
Bw=Bw+k*Q55(i)*I0(i);
Atheta=Atheta+Q11(i)*I1(i);
Btheta=Btheta-k*Q55(i)*I0(i);
C1theta=C1theta+Q11(i)*I2(i);
C2theta=C2theta-k*Q55(i)*I0(i);
Mu=Mu+Q11(i)*I1(i);
Mtheta=Mtheta+Q11(i)*I2(i);
end

% Du1=e31*I0(1)/esp(1);
% Du2=e31*I0(2)/esp(2);
% Dtheta1=e31*I1(1)/esp(1);
% Dtheta2=e31*I1(2)/esp(2);

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

S_total(1:n,1:n)=Au*D2AA;
S_total(1:n,n+1:2*n)=0;
S_total(1:n,2*n+1:3*n)=Cu*D2AA;

A_total(1:n,1:n)=J0*AA;
A_total(1:n,n+1:2*n)=0;
A_total(1:n,2*n+1:3*n)=J1*AA;

S_total(n+1:2*n,1:n)=0;                    
S_total(n+1:2*n,n+1:2*n)=Bw*D2AA;
S_total(n+1:2*n,2*n+1:3*n)=Bw*DAA;

A_total(n+1:2*n,1:n)=0;                    
A_total(n+1:2*n,n+1:2*n)=J0*AA;
A_total(n+1:2*n,2*n+1:3*n)=0;

S_total(2*n+1:3*n,1:n)=Atheta*D2AA;
S_total(2*n+1:3*n,n+1:2*n)=Btheta*DAA;
S_total(2*n+1:3*n,2*n+1:3*n)=C1theta*D2AA+C2theta*AA;

A_total(2*n+1:3*n,1:n)=J1*AA;
A_total(2*n+1:3*n,n+1:2*n)=0;
A_total(2*n+1:3*n,2*n+1:3*n)=J2*AA;

%CONDICOES DE FRONTEIRA

switch cfstr
    
    case {'ss'}
    
%equacao de fronteira - CASO SIMPLESMENTE APOIADO
b=find(x==0 | x==L );
% S_total(b,:)=[Au*DAA(b,:),zeros(2,n),Cu*DAA(b,:)]; %PARA NX
S_total(b,:)=[AA(b,:),zeros(2,n),zeros(2,n)];         %PARA U
S_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
S_total(b+2*n,:)=[Mu*DAA(b,:),zeros(2,n),Mtheta*DAA(b,:)]; %PARA M

A_total(b,:)=zeros(2,3*n);
A_total(b+n,:)=zeros(2,3*n);
A_total(b+2*n,:)=zeros(2,3*n);

    case {'cc'}        
%equacao de fronteira - CASO ENCASTRADO
b=find(x==0 | x==L );
S_total(b,:)=[AA(b,:),zeros(2,n),zeros(2,n)];          %PARA U
S_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
S_total(b+2*n,:)=[zeros(2,n),zeros(2,n),AA(b,:)];                 %PARA THETA

A_total(b,:)=zeros(2,3*n);
A_total(b+n,:)=zeros(2,3*n);
A_total(b+2*n,:)=zeros(2,3*n);

    case {'cl'}
b1=find(x==0); 
S_total(b1,:)=[AA(b1,:),zeros(1,n),zeros(1,n)];          %PARA U
S_total(b1+n,:)=[zeros(1,n),AA(b1,:),zeros(1,n)];       %PARA W
S_total(b1+2*n,:)=[zeros(1,n),zeros(1,n),AA(b1,:)];                 %PARA THETA

A_total(b1,:)=zeros(1,3*n);
A_total(b1+n,:)=zeros(1,3*n);
A_total(b1+2*n,:)=zeros(1,3*n);

b2=find(x==L);
S_total(b2,:)=[Au*DAA(b2,:),zeros(1,n),Cu*DAA(b2,:)]; %PARA NX
S_total(b2+n,:)=[zeros(1,n),Bw*DAA(b2,:),Bw*AA(b2,:)];       %PARA QX
S_total(b2+2*n,:)=[Mu*DAA(b2,:),zeros(1,n),Mtheta*DAA(b2,:)];       %PARA MX

A_total(b2,:)=zeros(1,3*n);
A_total(b2+n,:)=zeros(1,3*n);
A_total(b2+2*n,:)=zeros(1,3*n);
end

%%
%%

%%

%% DINAMICA

[lambda_vec,lambda]=eig(S_total,A_total);

lambda=diag(lambda,0); 
[lambda,indice]=sort(lambda);
lambda_vec=lambda_vec(:,(indice(:)));
lambda_norm=sqrt(lambda)*L^2*sqrt(rho*h/(Q11*I2));

% lambda=real(lambda);
%  lambda(lambda(:)==-inf | lambda(:)<=0)=NaN;
% 
% lambda=sqrt(lambda)*L^2*sqrt(rho*A/(E*I));%lambda=sqrt(lambda);




m=3; E=6e10; I=I2(1)+I2(2); A=h; G=2.3e10;
sol_exacta=(m*pi/L)^2*sqrt((E*I)/(rho*A))*sqrt(1-(((m*pi/L)^2*E*I)/(k*G*A+(m*pi/L)^2*E*I)));
% sol_exacta_norm=sol_exacta*L^2*sqrt(rho*A/(E*I));
%-------------graficos 


p=3;
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

figure(1)
plot(x,lambda_mode_w);
plot(x, lambda_mode_w(:,p));hold on;title(['Forma natural' num2str(p)]);legend(['wnorm(' num2str(p) ')=' num2str(lambda_norm(p),'%6.4f') ', w(' num2str(p) ')=' num2str(sqrt(lambda(p))/(2*pi),'%6.4f') 'Hz' ', w exata=' num2str(sol_exacta/(2*pi),'%6.4f') 'Hz' ,', erro=' num2str(100*(sol_exacta-sqrt(lambda(p)))/(sol_exacta),'%6.4f') '%'])


