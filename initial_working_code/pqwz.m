function [syse,np,nq,nw,nz] = pqwz(sys,ne,labeltype)

% [syse,np,nq,nw,nz] = pqwz(sys) creates extended system syse such that sys
% is the nominal part of the extended system syse that has inputs w and
% outputs z and extends the systems with "ne" extended inputs p1 ... pne and outputs q1 ... qne.
%
% Caution: ONLY WORKS FOR SYSTEMS WITH sys.D = zeros(no,ni)!
%
% General use of this kind of system is in robust control where w and z are
% performance channels and p1 ... pne and q1 ... qne are disturbance chanels on which the system uncertainty is connected.
%
% The resulting system syse is extended by means of truncating original
% iputs and outputs such that the syse.B = [sys.B, sys.B(1), ... sys.B(ne)] and syse.C = [sys.C; sys.C(1) ... sys.C(ne)].
%
% np, nq, nw and nz are number of inputs and outputs on each channel, p, q,
% w, z respectively.
%
% See 'help dynlbl" for type of labelign. Choose labeltype 1 or 2.

[nz, nw] = iosize(sys);
A = sys.a;
B = sys.b;
C = sys.c;
for j = 1:ne
    B = horzcat(B,sys.b);
    C = vertcat(C,sys.c);
end

D = zeros(size(C,1),size(B,2));
syse = ss(A,B,C,D);

if labeltype == 1
    syse.StateName = dynlbl('x',size(syse.A,1),labeltype);
    syse.InputName = vertcat(dynlbl('w',nw,labeltype),dynlbl('p',[ne,nw],labeltype+2)); % if you add new labeltypes then check this condition labeltype+2!!!
    syse.OutputName = vertcat(dynlbl('z',nz,labeltype),dynlbl('q',[ne,nz],labeltype+2));
elseif labeltype == 2
    syse.StateName = dynlbl('x',size(syse.A,1),labeltype);
    syse.InputName = vertcat(dynlbl('w',nw,labeltype),dynlbl('p',[ne,nw],labeltype+2));
    syse.OutputName = vertcat(dynlbl('z',nz,labeltype),dynlbl('q',[ne,nz],labeltype+2));
else
    error('labeltype should be 1 or 2')
end
[no, ni] = iosize(syse);
np = ni - nw;
nq = no - nz;

% check if maybe:
% syse.InputName = vertcat(dynlbl('p',[ne,nw],4),dynlbl('w',nw,2));
% syse.InputName = vertcat(dynlbl('q',[ne,nz],4),dynlbl('z',nz,2));
% this is how it is in Veenman paper Figure 2, lower right block is nominal...
