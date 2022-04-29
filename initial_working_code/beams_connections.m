if exist('beam_min_list','var') && ~isempty('beam_min_list')
    disp('Beams found. Continuing.')
else
    disp('No beam models detected. Creating beams...')
    zzz = 'zzz';
    eefsoib; % only used for comparison of the closed loops models obtained
    % with K matrix and connected model using connect -- they are equal
    %     beams;
    disp('Beams created.')
end
% ssred = beam_red;
G1 = beam_min_list(:,:,end); % finest discretization without reduction
% G1 = beam_min_list(1:2,1:2,end); % extracting subsystem
G1.InputName = 'w1';
G1.OutputName = 'z1';

bro1 = 6; % number of order of reduced model of system G1
G1t = ssred(:,:,bro1,end); % finest discretization with number of states bro1
% G1t = ssred(1:2,1:2,bro1,end); % extracting subsystem
G1t.InputName = 'w1';
G1t.OutputName = 'z1';


G2 = beam_min_list(:,:,end); % finest discretization without reduction
% G2 = beam_min_list(1:2,1:2,end); % extracting subsystem
G2.InputName = 'w2';
G2.OutputName = 'z2';

bro2 = 3; % number of order of reduced model of system G2
G2t = ssred(:,:,bro2,end);
% G2t = ssred(1:2,1:2,bro2,end); % extracting subsystem
G2t.InputName = 'w2';
G2t.OutputName = 'z2';

if exist('c1','var') && exist('c2','var') && exist('k1','var') && exist('k2','var')
    disp('Springs and dampers found. Continuing.')
else
    disp('No springs and dampers detected. Applying properties...')
    % Dampers' viscous coefficients
    c1 = 1e-2; % [Ns/m] [kg/s]
    c2 = 1e-1; %

    % Springs' stiffness
    k1 = 1e4; % [N/m] [kg/s^2]
    k2 = 1e5; %
    disp('Springs and dampers created.')
end

%% Springs-dampers dynamical system
sprdam = ss([-c1  0  -k1  0;
             0  -c2  0  -k2]);
% sprdam = ss([-c1  0 ;
%     0  -c2]); % extracting subsystem
sprdam.InputName = 's';
sprdam.OutputName = 'd';

% S2 = ss([-c1  0 ;
%     0  -c2]); % extracting subsystem
S2 = sprdam;
S2.InputName = 's2';
S2.OutputName = 'd2';

G3 = beam_min_list(:,:,end); % finest discretization without reduction
% G3 = beam_min_list(1:2,1:2,end); % extracting subsystem
G3.InputName = 'w3';
G3.OutputName = 'z3';

bro3 = 3; % number of order of reduced model of system G3
G3t = ssred(:,:,bro3,end);
% G3t = ssred(1:2,1:2,bro3,end); % extracting subsystem
G3t.InputName = 'w3';
G3t.OutputName = 'z3';

%%
% Additive model of difference between original and the reduced system G1
G1_G1t = G1 - G1t;
G1_G1t.InputName = 'q11';
G1_G1t.OutputName = 'p11';

% Additive model of difference between original and the reduced system G2
G2_G2t = G2 - G2t;
G2_G2t.InputName = 'q22';
G2_G2t.OutputName = 'p22';

% Additive model of difference between original and the reduced system G3
G3_G3t = G3 - G3t;
G3_G3t.InputName = 'q33';
G3_G3t.OutputName = 'p33';

%% Additive uncertainty filters
% Making input filter W11 that covers frequency response of local subsystem
w1 = logspace(1,5,200); % frequency range of interest
W11frd = frd(G1_G1t,w1); % vs reduced
W11frda = abs(W11frd); % jel mi to treba??
ord1 = 3; % order of filter W11 that approximates additive model G1 - G1t
[no1,ni1] = iosize(W11frda);
W11 = rss(ord1,no1,ni1); % preallocation of filter W1 approximating connected enviorement
for j = 1:no1
    for k = 1:ni1
        C1.UpperBound = [];
        C1.LowerBound = W11frda(j,k);
        b1 = fitmagfrd(W11frda(j,k),ord1,[],W11frda(j,k),C1);
        W11(j,k) = b1;
    end
end
W11.InputName = 'q11';
W11.OutputName = 'q1';

% Delta1 = ss(eye(ni1,no1));
Delta1 = ultidyn('Delta1',[ni1,no1]); % with uncertainty
Delta1.InputName = 'q1';
Delta1.OutputName = 'p1';

W12 = ss(eye(no1,ni1));
W12.InputName = 'p1';
W12.OutputName = 'p11';

% Making input filter that W21 covers frequency response of connected
% subsystem
w2 = logspace(1,5,200); % frequency range of interest
W21frd = frd(G2_G2t,w2); % vs reduced
W21frda = abs(W21frd); % jel mi to treba??
ord2 = 2; % order of filter W21 that approximates additive model G2 - G2t
[no2,ni2] = iosize(W21frda);
W21 = rss(ord2,no2,ni2); % preallocation of filter W1 approximating connected enviorement
for j = 1:no2
    for k = 1:ni2
        C2.UpperBound = [];
        C2.LowerBound = W21frda(j,k);
        b2 = fitmagfrd(W21frda(j,k),ord2,[],W21frda(j,k),C2);
        W21(j,k) = b2;
    end
end
W21.InputName = 'q22';
W21.OutputName = 'q2';

% Delta2 = ss(eye(ni2,no2));
Delta2 = ultidyn('Delta2',[ni2,no2]); % with uncertainty
Delta2.InputName = 'q2';
Delta2.OutputName = 'p2';

W22 = ss(eye(no2,ni2));
W22.InputName = 'p2';
W22.OutputName = 'p22';

%% Making connections
% summing junctions
s1 = sumblk('w1 = w + d',ni1);
% s2 = sumblk('z = z1 + p11',2); %4); % extracting subsystems
% s3 = sumblk('s = z + z22',2); %4); % extracting subsystems
% s4 = sumblk('z22 = p22 + z2',2); %4); % extracting subsystems
s2 = sumblk('z = z1 + p11',no1);
s3 = sumblk('s = z + z22',no1);
s4 = sumblk('z22 = p22 + z2',no1);
% forking points
f1 = sumblk('q11 = w1',ni1);
f2 = sumblk('w2 = d',ni1);
f3 = sumblk('q22 = w2',ni1);


%% Constructing closed loop systems

% Original system (G1 and G2 interconnected with beams and springs)
% creating connections (summing junctions) for original interconnected
% system
s1o = sumblk('w1 = w + d',ni1);
% s2o = sumblk('s = z + z2',2); %4); % extracting subsystem
% f1o = sumblk('z = z1',2); %4); % extracting subsystem
s2o = sumblk('s = z + z2',no1);
f1o = sumblk('z = z1',no1);
f2o = sumblk('w2 = d',ni1);
Gcl = connect(G1,G2,sprdam,s1o,s2o,f1o,f2o,'w','z');
% bodemag(Gcl,':k.',cl_ref,':r',{1e1 1e5}) % for testing, use with eefsoib

% Reduced system (G1t and G2t interconnected with beams and springs and weighted using additive uncertainty model)
Gclr = connect(G1t,W11,Delta1,W12,sprdam,G2t,W21,Delta2,W22,s1,s2,s3,s4,f1,f2,f3,'w','z');
% bodemag(Gclr,':k.',Gcl,':r',{1e1 1e5})

% Reducing original system to the same order of additive reduced order
balred_opts = balredOptions('StateProjection','truncate','ErrorBound','absolute');
[Gclbr,Gclbr_info] = balred(Gcl,size(Gclr.A,1),balred_opts); % series of reduced second subsystems
% bodemag(Gclr,':k.',Gclbr,':r',{1e1 1e5})
% bodemag(Gcl,':k.',Gclbr,':r',{1e1 1e5})

%% Creating weights for seconds system (using rest of the system)
s1e = sumblk('w1 = w + d',ni1);
s2e = sumblk('s = z + wx',no1);
f1e = sumblk('zx = d',ni1);
f2e = sumblk('z = z1',no1);
G1e = connect(G1,sprdam,s1e,s2e,f1e,f2e,{'wx','w'},{'zx','z'});
G11e = G1e('zx','w');
G11e_test = (sprdam*G1)*feedback(eye(2,2),sprdam*G1,+1);
% norm(G11e,Inf), norm(G11e_test,Inf)
% bodemag(G11e,'-k',G11e_test,'r.')
G12e = G1e('z','wx');
% G12e_test = (G1*sprdam)*feedback(eye(2,2),G1*sprdam,+1);
% norm(G12e,Inf), norm(G12e_test,Inf)
% bodemag(G12e,'-k',G12e_test,'r.')

% G1etest = connect(blksys_test,)
%% Creating weights for third system (using rest of the system)
s1e1 = sumblk('w1 = w + d',ni1);
s2e1 = sumblk('s = z1 + z2',no1);
s3e1 = sumblk('w2 = d + d2',ni1);
s4e1 = sumblk('s2 = z2 + wx2',no1);
f1e1 = sumblk('zx2 = d2',ni1);
f2e1 = sumblk('z = z1',no1);
G1e1 = connect(G1,sprdam,G2,S2,s1e1,s2e1,s3e1,s4e1,f1e1,f2e1,{'wx2','w'},{'zx2','z'});
G11e1 = G1e1('zx2','w');
% G11e1_test = ????;
% norm(G11e1,Inf), norm(G11e1_test,Inf)
% bodemag(G11e,'-k',G11e_test,'r.')
G12e1 = G1e1('z','wx2');
% G12e1_test = ????;
% norm(G12e1,Inf), norm(G12e1_test,Inf)
% bodemag(G12e1,'-k',G12e1_test,'r.')

% Making input filter W11e on the second system input that covers frequency of the first system and spring-dampers
w11e = logspace(1,5,200); % frequency range of interest
W11efrd = frd(G11e,w11e); % vs reduced
W11efrda = abs(W11efrd); % jel mi to treba??
ord11e = 2; % order of filter W21 that approximates additive model G2 - G2t
[no11e,ni11e] = iosize(W11efrda);
W11e = rss(ord11e,no11e,ni11e); % preallocation of filter W1 approximating connected enviorement
for j = 1:no11e
    for k = 1:ni11e
        C11e.UpperBound = [];
        C11e.LowerBound = W11efrda(j,k);
        b11e = fitmagfrd(W11efrda(j,k),ord11e,[],W11efrda(j,k),C11e);
        W11e(j,k) = b11e;
    end
end
% W11e.InputName = 'w';
% W11e.OutputName = 'w2';
% Reducing filter...
% [~,W11e_info] = balred(W11e,balred_opts);
% W11e_hsv = W11e_info.HSV; % Hankel singular values
% W11e_red_ord = sum(abs(W11e_hsv) > 10e-8);
% W11e_min = balred(W11e, W11e_red_ord, balred_opts); % minimal order beam state space model
% W11e_minimal = minreal(W11e_min); % checking if any more states can be discarded

% Making input filter W12e on the second system output that covers frequency of the first system and spring-dampers
w12e = logspace(1,5,200); % frequency range of interest
W12efrd = frd(G12e,w12e); % vs reduced
W12efrda = abs(W12efrd); % jel mi to treba??
ord12e = 2; % order of filter W21 that approximates additive model G2 - G2t
[no12e,ni12e] = iosize(W12efrda);
W12e = rss(ord12e,no12e,ni12e); % preallocation of filter W1 approximating connected enviorement
for j = 1:no12e
    for k = 1:ni12e
        C12e.UpperBound = [];
        C12e.LowerBound = W12efrda(j,k);
        b12e = fitmagfrd(W12efrda(j,k),ord12e,[],W12efrda(j,k),C12e);
        W12e(j,k) = b12e;
    end
end
% W11e.InputName = 'w';
% W11e.OutputName = 'w2';
% Reducing filter...
% [~,W12e_info] = balred(W12e,balred_opts);
% W12e_hsv = W12e_info.HSV; % Hankel singular values
% W12e_red_ord = sum(abs(W12e_hsv) > 10e-8);
% W12e_min = balred(W12e, W12e_red_ord, balred_opts); % minimal order beam state space model
% W12e_minimal = minreal(W12e_min); % checking if any more states can be discarded

G2te = W12e*G2t*W11e;
G2te.InputName = 'w2e';
G2te.OutputName = 'z2e';

G2e = W12e*G2*W11e;
G2e.InputName = 'w2e';
G2e.OutputName = 'z2e';

G2e_G2te = G2e - G2te;
% Making input filter that W21e covers frequency response of connected
% subsystem
w2e = logspace(1,5,200); % frequency range of interest
W21efrd = frd(G2e_G2te,w2e); % vs reduced
W21efrda = abs(W21efrd); % jel mi to treba??
ord2e = 2; % order of filter W21 that approximates additive model G2 - G2t
[no2e,ni2e] = iosize(W21efrda);
W21e = rss(ord2e,no2e,ni2e); % preallocation of filter W1 approximating connected enviorement
for j = 1:no2e
    for k = 1:ni2e
        C2e.UpperBound = [];
        C2e.LowerBound = W21efrda(j,k);
        b2e = fitmagfrd(W21efrda(j,k),ord2e,[],W21efrda(j,k),C2e);
        W21e(j,k) = b2e;
    end
end
W21e.InputName = 'q22';
W21e.OutputName = 'q2';

s1e = sumblk('z = z2e + p22',no12e);
f1e = sumblk('w2e = w',ni11e);
f2e = sumblk('q22 = w',ni11e);
Gclre = connect(G2te,W21e,Delta2,W22,s1e,f1e,f2e,'w','z'); %%%%% IS THIS AN ERROR?!?!??
%%%%% Ok so logic behind the upper statement is that I have a new system
%%%%% with external inputs 'w' and 'z' now acting on the second system, how
%%%%% does this play out in the grand scheme? when I plug in newly obtained
%%%%% (refined) weighting filters into aggregate model? Ok so I check it
%%%%% later with Gclre4, so this should be fine
% test with the same filter W21 as before but with weighted G2e
Gclre1 = connect(G2te,W21,Delta2,W22,s1e,f1e,f2e,'w','z');
% bodemag(Gcl,'k',Gclr.NominalValue,'--r',Gclre.NominalValue,':b',Gclre1.NominalValue,':m')

% test with the W21s static filters
nW21s1 = norm(W21e,Inf);
W21s1 = ss([nW21s1, nW21s1; nW21s1, nW21s1]);
W21s1.InputName = 'q22';
W21s1.OutputName = 'q2';
Gclre2 = connect(G2te,W21s1,Delta2,W22,s1e,f1e,f2e,'w','z');
% test with the W21s static filters
nW21s2 = 0.01; % approximate lower limit of W21
W21s2 = ss([nW21s2, nW21s2; nW21s2, nW21s2]);
W21s2.InputName = 'q22';
W21s2.OutputName = 'q2';
Gclre3 = connect(G2te,W21s2,Delta2,W22,s1e,f1e,f2e,'w','z');

% bodemag(W21e,'k',W21,'r',W21s1,'--g',W21s2,'--r')

% Reduced system (G1t and G2t interconnected with beams and springs and weighted using additive uncertainty model)
% for checking the newly obtained refined weights
Gclre4 = connect(G1t,W11,Delta1,W12,sprdam,G2t,W21s1,Delta2,W22,s1,s2,s3,s4,f1,f2,f3,'w','z');


disp('Interconnected beams created.')
%% IQC analysis

% Generating open loop nominal plant
if class(Gclre) == char('uss')
    Me = lftdata(Gclre);
    Me1 = lftdata(Gclre1);
    Me2 = lftdata(Gclre2);
    Me3 = lftdata(Gclre3);
    Me4 = lftdata(Gclre4);
    if (sum(real(eig(Me)) >= 0))
        error('Nominal open loop plant Me is not stable.')
    end
else
    error('Nominal open loop plant Me can only be obtained from uncertain state-space models using lftdata.')
end

if class(Gclr) == char('uss')
    M = lftdata(Gclr);
    if (sum(real(eig(M)) >= 0))
        error('Nominal open loop plant M is not stable.')
    end
else
    error('Nominal open loop plant M can only be obtained from uncertain state-space models using lftdata.')
end
disp('Stable open loop plants created.')

return;

clear ga de1 de2 de pe prob prob1 prob2 prob3 prob4
%% IQC analysis
abc = 1;
if abc == 1 % original IQC analysis
    % % Define uncertainty block Delta1
        de1 = iqcdelta('de1','InputChannel',[1:2],'OutputChannel',[1:4],'StaticDynamic','D','NormBounds',0.02,'Structure','FB');
%     de1 = iqcdelta('de1','InputChannel',[1:2],'OutputChannel',[1:2],'StaticDynamic','D','NormBounds',0.02,'Structure','FB');
    % % Assign IQC-multiplier to uncertainty block Delta1
    % de1 = iqcassign(de1,'ultid','Length',3,'PoleLocation',-1);


    % % Define uncertainty block Delta2
        de2 = iqcdelta('de2','InputChannel',[3:4],'OutputChannel',[5:8],'StaticDynamic','D','NormBounds',0.02,'Structure','FB');
%     de2 = iqcdelta('de2','InputChannel',[3:4],'OutputChannel',[3:4],'StaticDynamic','D','NormBounds',0.02,'Structure','FB');
    % % Assign IQC-multiplier to uncertainty block Delta2
    % de2 = iqcassign(de2,'ultid','Length',3,'PoleLocation',-1);

    de = blkdiag('de',de1,de2);
    de = iqcassign(de,'ultid','Length',3,'PoleLocation',-1);

    % Define performance block (with two uncertainties)
        pe = iqcdelta('pe','ChannelClass','P','InputChannel',[5:6],'OutputChannel',[9:12],'PerfMetric','L2');
%     pe = iqcdelta('pe','ChannelClass','P','InputChannel',[5:6],'OutputChannel',[5:6],'PerfMetric','L2');


    % Define uncertainty block Delta - only second system with uncertainty
    % de = iqcdelta('de','InputChannel',1:2,'OutputChannel',1:2,'StaticDynamic','D','Structure','FB','NormBounds',1);
    % Assign IQC-multiplier to uncertainty block - only second system with uncertainty
    % de = iqcassign(de,'ultid','Length',3,'PoleLocation',-1);
    % Define performance block (with one uncertainty)
    % pe = iqcdelta('pe','ChannelClass','P','InputChannel',3:4,'OutputChannel',3:4,'PerfMetric','L2');

    tic
    % Perform IQC-analysis
    prob = iqcanalysis(M,{de,pe},'SolChk','on','eps',1e-8);
    toc

elseif abc == 2 % refinement of the uncertainty of the second system
    % % Define uncertainty block Delta1
    de = iqcdelta('de','InputChannel',[1:2],'OutputChannel',[1:2],'StaticDynamic','D','NormBounds',0.02,'Structure','FB');
    de = iqcassign(de,'ultid','Length',3,'PoleLocation',-1);

    % Define performance block (with two uncertainties)
    pe = iqcdelta('pe','ChannelClass','P','InputChannel',[3:4],'OutputChannel',[3:4],'PerfMetric','L2');

    tic
    % Perform IQC-analysis
    prob = iqcanalysis(Me,{de,pe},'SolChk','on','eps',1e-8);
    %     prob = iqcanalysis(Me,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
    toc

    tic
    prob1 = iqcanalysis(Me1,{de,pe},'SolChk','on','eps',1e-8);
    %     prob1 = iqcanalysis(Me1,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
    toc

    tic
    prob2 = iqcanalysis(Me2,{de,pe},'SolChk','on','eps',1e-8);
    %     prob2 = iqcanalysis(Me2,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
    toc

    tic
    prob3 = iqcanalysis(Me3,{de,pe},'SolChk','on','eps',1e-8);
    %     prob3 = iqcanalysis(Me3,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
    toc
elseif abc == 3 % checking if the newly obtained weigths give solution to the original problem
    % % Define uncertainty block Delta1
    de = iqcdelta('de','InputChannel',[1:4],'OutputChannel',[1:4],'StaticDynamic','D','NormBounds',1,'Structure','FB');
    de = iqcassign(de,'ultid','Length',1,'PoleLocation',-1);

    % Define performance block (with two uncertainties)
    pe = iqcdelta('pe','ChannelClass','P','InputChannel',[5:6],'OutputChannel',[5:6],'PerfMetric','L2');


    tic
    prob4 = iqcanalysis(Me4,{de,pe},'SolChk','on','eps',1e-8);
%         prob4 = iqcanalysis(Me4,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','sdpt3');
    toc
end

% bodemag(W21e,':k',W21,':r')

%% Test after prob
% norm(Me1('z','w'),Inf)