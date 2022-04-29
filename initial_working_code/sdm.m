function [kk, hh, rr,ck] = sdm(n,ni,no,c1,c2,k1,k2,type)
% syntax: [kk, hh, rr] = sdm(n,ni,no,c1,c2,k1,k2)
% sdm (short for spring-damper-matrix) is used for creating connection matrix kk for given number of beams
% which are mutually interconnected at two places with springs and dampers.
% There are two external forces acting on the first beam and thefirst and the
% last beam has springs and dampers only on one side.
% Input parametars are:
%                       n - number of beams, n >= 2 (i.e at least two beams
%                           need to be interconnected)
%                       ni - number of beam at which the external inputs are
%                           imposed (ni >= 1, i.e can be at the first beam of the
%                           interconnected system)
%                       no - number of beam whose outputs are used as
%                           external outputs (no => 1, i.e internal output of
%                           any of the beams, including the first one, can be
%                           used as external output)
%                       c1, c2 - viscous coefficients
%                       k1, k2 - spring stiffnesses
m = 2; % number of external inputs of interconnected dynamical system (forces) i.e [F1(t) F2(t)].'
p = 4; % number of external outputs of interconnected dynamical system (velocities of first
% and second node and displacaments of first and second node)
% i.e [v1(t) v2(t) x1(t) x2(t)].'
mj = 2; % number of internal inputs of dynamical system (forces on one beam) i.e [F1(t) F2(t)].'
pj = 4; % number of external outputs of the dynamical system (velocities of first
% and second node and displacaments of first and second node)
% i.e [v1(t) v2(t) x1(t) x2(t)].'

ns = n + n -1; % number of all systems, including spr-dam systems
nsd = ns - n; % number of spr-dam systems
msd = 4; % number of inputs to spr-dam systems
psd = 2; % number of outputs of spr-dam systems

switch type
    case 'type1'
        if n <= 1 || length(n) > 1
            error('Number of beams must be integer n >= 2')
        else
            ck = [-c1, 0, -k1, 0  ;
                0, -c2, 0, -k2];
            ckk = [ck, -ck;
                -ck,  ck];
            kk = zeros(n*mj,n*pj);
            for j = 1:n-1
                kk(mj * (j-1) + 1 : mj * (j+1), pj * (j-1) + 1: pj * (j+1)) = kk(mj * (j-1) + 1 : mj * (j+1), pj * (j-1) + 1: pj * (j+1)) + ckk; % add overlapping block diagonal matrices
            end
        end

        hh = zeros(n*mj,m); % original
        if ni <= 0 || ni > n
            error('0 < ni <= n')
        else
            hh(mj*(ni-1)+1:ni*mj,1:m) = eye(mj,m); % original
        end

        rr = zeros(p,n*pj); % original
        if no <= 0 || no > n
            error('0 < no <= n')
        else
            rr(1:pj,pj*(no-1)+1:no*pj) = eye(p,pj); % original
        end

    case 'type2'
        ck = ss([-c1, 0, -k1, 0  ;
            0, -c2, 0, -k2]);

        if n <= 1
            error('Number of beams must be n >= 2')
        else

            kk = [];
            for l = 1:ns % rows
                kkr = [];
                for k = 1:ns % columns
                    if mod(l,2) == 1 && mod(k,2) == 1 % beam-beam detected
                        if k == l + 1 || l == k + 1
                            kkr = horzcat(kkr,eye(mj,pj));
                        else
                            kkr = horzcat(kkr,zeros(mj,pj));
                        end
                    elseif mod(l,2) == 1 && mod(k,2) == 0
                        if k == l + 1 || l == k + 1
                            kkr = horzcat(kkr,eye(mj,psd));
                        else
                            kkr = horzcat(kkr,zeros(mj,psd));
                        end
                    elseif mod(l,2) == 0 && mod(k,2) == 1
                        if k == l + 1 || l == k + 1
                            kkr = horzcat(kkr,eye(msd,pj));
                        else
                            kkr = horzcat(kkr,zeros(msd,pj));
                        end
                    elseif mod(l,2) == 0 && mod(k,2) == 0
                        if k == l + 1 || l == k + 1
                            kkr = horzcat(kkr,eye(msd,psd));
                        else
                            kkr = horzcat(kkr,zeros(msd,psd));
                        end
                    end
                end
                kk = vertcat(kk,kkr);
            end

        end

        hh = zeros(n*mj+nsd*msd,m);
        if ni <= 0 || ni > n
            error('0 < ni <= n')
        else
            hh((mj+msd)*(ni-1)+1:(ni-1)*(mj+msd)+1+mj-1,1:m) = eye(mj,m);
        end

        rr = zeros(p,n*pj+nsd*psd);
        if no <= 0 || no > n
            error('0 < no <= n')
        else
            rr(1:p,(pj+psd)*(no-1)+1:(no-1)*(pj+psd)+1+pj-1) = eye(p,pj);
        end
end



