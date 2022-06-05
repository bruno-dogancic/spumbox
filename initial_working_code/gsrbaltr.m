function [sysr,W,T,nr,St] = gsrbaltr(sys,ord)

% Algorithm 3.1 (ASOMROCS) Generalized square root balanced truncation method
% syntax: [sysr,W,T,nr,St] = gsrbaltr(sys,ord)
% 1. Compute the lower Cholesky factors of contr. and obs. grammians
% 2. Compute the SVD
% 3. Compute the reduced system
% Outputs are: sysr - Reduced system of order nr = (order(sys) - ord)
%              W, T - diagnonal projection matrices
%              St - truncated singular values

P = gram(sys,'c');
Q = gram(sys,'o');
Lp = chol(P,'lower');
Lq = chol(Q,'lower');
[U,S,V] = svd(Lp'*(sys.e)'*Lq);
U1 = U(:,1:ord);
S1 = S(1:ord,1:ord);
St = diag(S(ord+1:end,ord+1:end));
V1 = V(:,1:ord);
W = Lq*V1/sqrt(S1);
T = Lp*U1/sqrt(S1);
Et = W'*sys.e*T;
At = W'*sys.a*T;
Bt = W'*sys.b;
Ct = sys.c*T;
rc = size(sys);
sysr = dss(At,Bt,Ct,zeros(rc),Et);
nr = size(At,1);
