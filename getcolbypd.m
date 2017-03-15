function [colset,lam]=getcolbypd(mat,rhs)
%  
% Apply adaptive greedy to solve consistent but highly redundant RBF system 
%     (mat)*(lam)=(rhs) 
% for approximated solution (lam)
% 
% Algorithm (without matrix-free feature):
% L. Ling and R. Schaback. An improved subspace selection algorithm for 
% meshless collocation methods. International Journal for Numerical Methods
% in Engineering. 
% http://dx.doi.org/10.1002/nme.2674
% 
%     input: MxN (mat) with M>=N
%     input: Mx1 (rhs) 
%     output: colset = indices of columns in (mat) selected
%     output: lam = Nx1 solution to the system (optional)
% 
% To obtain reduced system (A_r)*(lam_r)=rhs to the original,
% >> A_r=A(:,colset);
% >> lam_r=Ar\rhs;
% 
% Test problem:
% >> n=500;
% >> A=rand(n)^10; % make A ill-conditioned
% >> b=rand(n,1);  
% >> norm(A*(A\b)-b,'inf') % check residual error by \
% >> [colset,lam]=getcolbypd(A,b); 
% 
% by Leevan Ling (lling@hkbu.edu.hk), 2009.

[m n]=size(mat);
[mb nb]=size(rhs);
if nb~=1, rhs=rhs'; [mb nb]=size(rhs); end
if mb~=m, error('Input dimensions do not match'); end
[mv imv]=max(abs(rhs));        % get maximum residual
[nv inv]=max(abs(mat(imv,:))); % take maximum in correspondent row
x=zeros(n,1);
x(inv,1)=rhs(imv,1)/mat(imv,inv);
colset=[inv]; % this will be the selected column set
rowset=[imv]; % this will be the selected row set
v=zeros(m,1);
v(imv,1)=-x(inv,1)/mat(imv,inv);
% now loop for the primal-dual algorithm
% no inverse update implemented
lastwarn('');
for k=2:min(m,n)
    r=mat*x-rhs;
    q=x+mat'*v;
    [mv imv]=max(abs(r));
    [nv inv]=max(abs(q));
    rowset=[rowset imv];
    colset=[colset inv];    
    matloc=mat(rowset,colset);
    [lastmsg, lastid] = lastwarn;
    if ~isempty(lastmsg), break; end
    rhsloc=rhs(rowset);
    xloc=matloc\rhsloc;
    vloc=-(matloc')\xloc;
    x=zeros(n,1);
    x(colset,1)=xloc;
    v=zeros(m,1);
    v(rowset,1)=vloc;
end
colset=unique(colset);
disp(sprintf('Adaptive greedy algorithm stopped after selecting %g columns',length(colset)))
if nargout>1
    lam = zeros(n,1);
    warning off
    lam(colset) = mat(:,colset)\rhs;
    warning on
    disp(sprintf('inf-norm of residual = %g, 2-norm of residual = %g',...
        norm(mat(:,colset)*lam(colset)-rhs,'inf'),...
        norm(mat(:,colset)*lam(colset)-rhs)))
end


