function [nags,Ls,Xs,sol]=nh_fun(r0,rN,W,options)
% SDP computation of Nagaoka--Hayashi CR-bound for multi-parameter estimation [1].
%
% [nags, Ls, Xs, sol ] = nag_fun(S0,{S1,S2,S3},W,opts)
%
% S0 is the state and {Sj} are the state derivatives
%
% W: third parameter (optional) is the weight matrix
%
% fourth parameter (optional) 
%   opt=sdpsettings('verbose',1,'debug',1,'solver','mosek');
%
% references:
% [1] Conlon, L. O., Suzuki, J., Lam, P. K., & Assad, S. M. (2021).
% Efficient computation of the Nagaoka–Hayashi bound for multiparameter
% estimation with separable measurements. npj Quantum Information, 7(1)
% arXiv:2008.02612

dim=size(rN{1},1);
n=length(rN);

if nargin==2
    W=eye(n);
    options=sdpsettings('verbose',1,'debug',1,'solver','sedumi');
end

if nargin==3
    options=sdpsettings('verbose',1,'debug',1,'solver','sedumi');
end

% create hermitian submatrices L_jk
for j=1:n
    v{j,j}= sdpvar(dim,dim,'hermitian','complex');
    for k= j+1:n
        v{j,k}= sdpvar(dim,dim,'hermitian','complex');
    end
end

% assemble submatrices into a symmetric L matrix.
v8=[];
for j=1:n
    v8row=[];
    for k=1:n
        v8row=[v8row,v{min(j,k),max(j,k)}];
    end
    v8=[v8;v8row];
end

% create hermitian submatrices X_j
for j=1:n
    Xcell{j}= sdpvar(dim,dim,'hermitian','complex');
end

% assemble submatrices into X.
X=[];Xt=[];
for j=1:n
    X=[X;Xcell{j}];
    Xt=[Xt,Xcell{j}];
end

Y=[v8,X;Xt,eye(dim)];
constraint1=[Y>=0];

constraint2a=[];
constraint2b=[];
for j=1:n
    constraint2a=[constraint2a,trace(r0*Xcell{j})==0];
    for k=1:n
        constraint2b=[constraint2b,trace(rN{k}*Xcell{j})==(j==k)];
    end
end
constraints=[constraint1,constraint2a,constraint2b];

rBold=kron(W,r0);
objective= trace(rBold*v8);

sol=optimize(constraints,objective,options);

Ls=value(v8);
Xs=value(X);
nags=value(objective);
