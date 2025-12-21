

s=1:100;
vv=zeros(size(s));
rhoBSOR=zeros(size(s));
rhoBJ=zeros(size(s));
rhoBGS=zeros(size(s));
indice=0;
a=-2;
b=2;
t=linspace(a,b,100);

n=3;

for i=1:100
    indice=indice+1;
    vv(indice)=t(i);
    A=es336_punto1(t(i));
    D=diag(diag(A));
% costruzione della matrice L
L=tril(A,-1);
% costruzione della matrice U
U = triu(A,1);
% costruzione della matrice B1
omega=0.1; %%Ho scelto questo valore di omega così per Ostrowski_Reich converge.
BSOR = (D+omega*L)\((1-omega)*D-omega*U);
z0=ones(n,1);

    rhoBSOR(indice)=max(abs(eigs(BSOR)));

    BGS =-(D+L)\U;

    rhoBGS(indice)=max(abs(eigs(BGS)));
    invD=diag(1./diag(A));  
    
    BJ=-invD*(L+U); 

    rhoBJ(indice)=abs(potenze(BJ,z0,1000,1e-6,2));







end 

figure(1)
semilogy(vv,rhoBJ,"Red",vv,rhoBGS,"Blue",vv,rhoBSOR,"Black")
legend("rhoBJ","rhoBGS","rhoBSOR")
title("Confronto raggi spettrali")
xlabel("t")
ylabel("rho")
grid on 