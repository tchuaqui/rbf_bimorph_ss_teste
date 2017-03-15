function [c_opt] = cross_validation(cfstr,n,L,x,Au,Cu,Bw,Atheta,Btheta,C1theta,C2theta,vetor_carga,Du1,Du2,Dtheta1,Dtheta2,phi1,phi2,Mu,Mtheta)
contador=1;
for c=0.01:0.1:10   
[xi,xj]=meshgrid(x,x);
AA=g(c,xi,xj);
DAA=dgdx(c,xi,xj);
D2AA=d2gdx2(c,xi,xj);
% 
% aphi1=(AA^-1)*phi1;
% aphi2=(AA^-1)*phi2;

S_total(1:n,1:n)=Au*D2AA;
S_total(1:n,n+1:2*n)=0;
S_total(1:n,2*n+1:3*n)=Cu*D2AA;

S_ele(1:n,1:n)=Du1*DAA;
S_ele(1:n,n+1:2*n)=Du2*DAA;

S_total(n+1:2*n,1:n)=0;                    
S_total(n+1:2*n,n+1:2*n)=Bw*D2AA;
S_total(n+1:2*n,2*n+1:3*n)=Bw*DAA;

S_ele(n+1:2*n,1:n)=0;
S_ele(n+1:2*n,n+1:2*n)=0;
 
S_total(2*n+1:3*n,1:n)=Atheta*D2AA;
S_total(2*n+1:3*n,n+1:2*n)=Btheta*DAA;
S_total(2*n+1:3*n,2*n+1:3*n)=C1theta*D2AA+C2theta*AA;

S_ele(2*n+1:3*n,1:n)=Dtheta1*DAA;
S_ele(2*n+1:3*n,n+1:2*n)=Dtheta2*DAA;
switch cfstr
    
    case {'ss'}
    
%equacao de fronteira - CASO SIMPLESMENTE APOIADO
b=find(x==0 | x==L );
% S_total(b,:)=[Au*DAA(b,:),zeros(2,n),Cu*DAA(b,:)]; %PARA NX
S_total(b,:)=[AA(b,:),zeros(2,n),zeros(2,n)];
S_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
S_total(b+2*n,:)=[Mu*DAA(b,:),zeros(2,n),Mtheta*DAA(b,:)]; %PARA M


    case {'cc'}        
%equacao de fronteira - CASO ENCASTRADO
b=find(x==0 | x==L );
S_total(b,:)=[AA(b,:),zeros(2,2*n)];          %PARA U
S_total(b+n,:)=[zeros(2,n),AA(b,:),zeros(2,n)];       %PARA W
S_total(b+2*n,:)=[zeros(2,2*n),AA(b,:)];                 %PARA THETA


    case {'cl'}
b1=find(x==0); 
S_total(b1,:)=[AA(b1,:),zeros(1,n),zeros(1,n)];          %PARA U
S_total(b1+n,:)=[zeros(1,n),AA(b1,:),zeros(1,n)];       %PARA W
S_total(b1+2*n,:)=[zeros(1,n),zeros(1,n),AA(b1,:)];                 %PARA THETA

b2=find(x==L);
S_total(b2,:)=[Au*DAA(b2,:),zeros(1,n),Cu*DAA(b2,:)]; %PARA NX
S_total(b2+n,:)=[zeros(1,n),Bw*DAA(b2,:),Bw*AA(b2,:)];       %PARA QX
S_total(b2+2*n,:)=[Mu*DAA(b2,:),zeros(1,n),Mtheta*DAA(b2,:)];       %PARA MX
end

e=repmat(eye(numel(x),numel(x)),3,1);  %CONSIDERANDO U, W E THETA
% e=eye(numel(x),numel(x));     %CONSIDERANDO APENAS W
for i=1:size(x,2) 
    m(:,i)=S_total(1:3*numel(x),1:3*numel(x))\e(:,i);   %CONSIDERANDO U, W E THETA
%     m(:,i)=S_total(n+1:2*n,n+1:2*n)\e(:,i);     %CONSIDERANDO W
end

vetor_feletrica=S_ele*[phi1;phi2];
vetor_feletrica(1)=0; vetor_feletrica(n)=0; vetor_feletrica(n+1)=0; 
vetor_feletrica(2*n)=0; vetor_feletrica(2*n+1)=0; vetor_feletrica(3*n)=0;
% switch cfstr
%     case {'ss'}
%       aux=Dtheta1*AA*phi1+Dtheta2*AA*phi2;
%      vetor_feletrica(b+2*n)=aux(b);
%     case {'cl'}
%       aux1=Du1*AA*phi1+Du2*AA*phi2;
%       aux2=Dtheta1*AA*phi1+Dtheta2*AA*phi2;
%       vetor_feletrica(b2)=aux1(b2);
%       vetor_feletrica(b2+2*n)=aux2(b2);
% end   


a=S_total\(vetor_carga-vetor_feletrica);
% a=a(n+1:2*n);
% for k=1:size(x,2) 
for k=1:size(x,2)
erro_rippa(k)=a(k)/m(k,k);
end
EE_rippa(contador)=norm(erro_rippa,2);
cc(contador)=c;
contador=contador+1;
end

[mx,indx]=min(EE_rippa);
c_opt=cc(indx);

end

