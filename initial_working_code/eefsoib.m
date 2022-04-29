%% Error estimates for series of interconnected beams
% Main ideas is to estimate an error introduced to the system by choosing coarse discretization
% and reducing the system.
% The system consists of a series of beams mutually interconnected with
% springs and dampers. The first beam in the series is considered correct
% (reference, exact or nominal model) and some of the beams that are
% interconnected are modeled using coarser finite element mesh but as well
% reduced in order to calculate error on the resulting system.
% Variables worth considering:
%   ns = number of interconnected systems (2 or more)
%   x = vector of integer values defining discretization, first one is the finest
%       discretization (i.e nd) since it represents the reference model
%   y = vector of integer values defining the order of reduction of the
%       system, first one is unreduced (i.e Inf) which corresponds to the reference
%       model, values can be from 1,2,3....,maxorder, Inf, where Inf is used to
%       define unreduced system

nb = 2;% % number of beams
ns = nb + nb -1; % number of systems (beams and springs between)
nd = 5;% % number of discretizations
x = [nd, nd-1];%, 10, 10];%, 10, 10, 10, 10, 10, 10];%
% x = [nd, nd];%, 10, 10, 10, 10, 10, 10, 10, 10];%
msord = 20; % maximum order of the reduced subsystems (if it cannot complete check hsv_tol below)
y = [Inf, msord-15];%, msord-14, msord-12];%, msord-10, msord-10, msord, msord, msord-7, msord-13];%
% y = [Inf, Inf];%, msord-14, msord-12, msord-10];%, msord-10, msord, msord, msord-7, msord-13];%

%% Cheking inputs
if nb < 2
    error('Number of beams nb >= 2.')
end
if (length(x) ~= nb || length(y) ~= nb)
    error ('Length of x and y must be nb = %i or since length(x) = length(y) = %i, nb must be nb = %i.',nb,length(x), length(y))
elseif (length(x) ~= length(y))
    error ('Length of x and y must be same')
end

%% Beam properties
e_beam = 210e9; % [GPa] [GN/m^2] [Gkg/s^2/m]
l_beam = 1; % [m]
d_beam = 0.01; % [m]
rho_beam = 7800; % [kg/m^3]

m_beam = rho_beam * l_beam * d_beam^2 * pi / 4;
i_beam = pi * d_beam^4 / 64;
ei_beam = e_beam * i_beam;

%% Series of discretized beams
nFE = zeros(1,nd); % number of diskretized models is ndm
for j = 1:length(nFE)
    nFE(j) = 9*j; % number of finite elements for each model 9, 18,..., 90
end

% Alocation for each variable that chages size in for loops in this section
% beam_ss_list(:,:,length(nke)) = rss(500,4,2);
% beam_dss_list(:,:,length(nke)) = beam_ss_list(:,:,length(nke));

% Creating series of finer discretizations
for j = 1:length(nFE)

    % discretized beam (mass and stiffness matrices)
    alpha = 0.02; % Coefficient of dynamic moment of inertia (rotation) for mass matrix
    [mat_m, mat_k] = be_beam_fe(nFE(j), m_beam, ei_beam, l_beam, alpha);

    % eigenvectors and eigenfrequencies for modal damping
    [eigvekt, eigfreq] = eig(mat_k, mat_m,'vector');
    omega = sort(sqrt(eigfreq));

    % damping matrix
    zeta = 0.05; % modal damping coefficient
    n_omega = 8;
    [ alfa, beta ] = rayparam(omega(1:n_omega), zeta);

    mat_p = alfa * mat_m + beta * mat_k;

    % Defining the input and output matrices
    nke3i = (nFE(j) / 3);

    % Input vector is [F1(t); F2(t)] - forces in nodes at l/3 and 2l/3
    mat_b1 = zeros(2 * nFE(j), 2);
    mat_b1(2 * nke3i, 1) = 1;
    mat_b1(4 * nke3i, 2) = 1;

    % Output vectors is [v1(t); v2(t); x1(t); x2(t)] - velocities and displacements
    mat_c1 = zeros([2, 2 * nFE(j)]);
    mat_c1(1, 2 * nke3i) = 1;
    mat_c1(2, 4 * nke3i) = 1;

    % Creating descriptor state-space dynamical system
    [mat_a_dss, mat_b_dss, mat_c_dss, mat_d_dss, mat_e_dss] = mehdss(mat_m, mat_p, mat_k, mat_b1, mat_c1);
    beam_dss_list(:,:,j) = dss(mat_a_dss,mat_b_dss,mat_c_dss,mat_d_dss,mat_e_dss);

    % Creating state-space dynamical system
    [mat_a, mat_b, mat_c, mat_d] = mehss(mat_m, mat_p, mat_k, mat_b1, mat_c1);
    beam_ss_list(:,:,j) = ss(mat_a, mat_b, mat_c, mat_d);
    %     beam_ss_list(:,:,j) = ss(beam_dss_list(:,:,j),'explicit'); %  alternative way to calculate state-space
end

%% Minimizing models
% Alocation for each variable that chages size in for loops in this section
% beam_ss_min_list(:,:,length(nke)) = greda_ss_list(:,:,length(nke));

% Creating treshold for truncating small hankel singular values in order to
% get as close to minimal realization as posible, but to leave some space
% for testing the algorithm
hsv_tol = 1e-12; % treshold for truncation
balred_opts = balredOptions('StateProjection','truncate','ErrorBound','absolute');

% Creating series of initially ("mildly") reduced discretized models
% which will be used as minimal reference model afterwards
for j = 1:length(nFE)
    [~,beam_red_info] = balred(beam_ss_list(:,:,j),balred_opts);
    beam_hsv = beam_red_info.HSV; % Hankel singular values
    red_ord = sum(abs(beam_hsv) > hsv_tol);
    beam_min_list(:,:,j) = balred(beam_ss_list(:,:,j), red_ord, balred_opts); % minimal order beam state space model
    beam_minimal_list(:,:,j) = minreal(beam_min_list(:,:,j)); % checking if any more states can be discarded
end

%% Creating reduced models of subsystems

for j = 1:length(nFE)
    maxrst = msord + 1; % max order of the system
    ssredord = 1:maxrst-1; % series of orders for model reduction of second system
    [ssred(:,:,:,j),bal_info(j)] = balred(beam_min_list(:,:,j),ssredord,balred_opts); % series of reduced second subsystems

    % H_inf norm of subsystem error, upper and lower bounds
    for jk = 1:length(ssredord)
        ssrederr(:,:,jk,j) = beam_min_list(:,:,j) - ssred(:,:,jk,j);
        ssrederrnorm(jk,j) = norm(ssrederr(:,:,jk,j),"inf");
        gammau(jk,j) = 2 * sum(bal_info(j).HSV(jk+1:end));
        %         gammau_test(jk,j) = sum(bal_info(j).HSV(jk+1:end));
        gammal(jk,j) = bal_info(j).HSV(jk+1); % for plotting only, no sense since first order reduced system doesn't have lower bound...
    end
end

% if isempty(zzz) == 0
%     return
% end
%% Spring and dampers parameters
% case 1
% viscous dampers
c1 = 1e-2; % [Ns/m] [kg/s]
c2 = 1e-6; %
% springs
k1 = 1e-8; % [N/m] [kg/s^2]
k2 = 1e-5; %

% % case 2
% viscous dampers
% c1 = 1e-1; % [Ns/m] [kg/s]
% c2 = 1e-1;
% % springs
% k1 = 1e+2; % [N/m] [kg/s^2]
% k2 = 1e+2;

%% Creating interconnected (closed-loop) model of ns beams in series coupled with springs and viscous dampers
% First subsystem is the finest discretization minimal system
% Other subsystems will be chosen with respect to the x and y

% Creating connection matrices 'type1'
ni = 1; % number of the system at which the external inputs are imposed
no = 1; % number of the system from which the outputs are also considered external outputs
[kk, hh, rr] = sdm(nb,ni,no,c1,c2,k1,k2,'type1');

pj = 4; % number of outputs of the beam subsystems
mj = 2; % number of inputs of the beam subsystems
pil = []; % PI_l matrix defining what systems are reduced (in terms of their outputs)
pir = []; % PI_r matrix defining what systems are reduced (in terms of their inputs)
ssdiag = []; % the diagonal system, first should be unreduced
ssdiag_ref = []; % diagonal reference system

for j = 1:nb
    ssdiag_ref = append(ssdiag_ref,beam_min_list(:,:,length(nFE)));
    %     ssdiag_ref = append(ssdiag_ref,beam_min_list(:,:,length(x(j)))); % this one also makes sense...
    if y(j) == Inf
        pil = blkdiag(pil,zeros(pj));
        pir = blkdiag(pir,zeros(mj));
        ssdiag = append(ssdiag,beam_min_list(:,:,x(j)));
    else
        pil = blkdiag(pil,eye(pj));
        pir = blkdiag(pir,eye(mj));
        ssdiag = append(ssdiag,ssred(:,:,y(j),x(j)));
    end
end

% Creating connection matrices 'type2'
[K, H, R, S] = sdm(nb,ni,no,c1,c2,k1,k2,'type2');
psd = 2; % number of outputs of the spring-damper subsystems
msd = 4; % number of inputs of the spring-damper subsystems
PIl = []; % PI_l matrix defining what systems are reduced (in terms of their outputs)
PIr = []; % PI_r matrix defining what systems are reduced (in terms of their inputs)
Ga = []; % the diagonal system, first should be unreduced
Ga_ref = []; % diagonal reference system

% Creating extended subsystem choosing vectors
ye = zeros(1,length(y)+length(y)-1);
xe = zeros(1,length(x)+length(x)-1);
for j = 1:nb
    ye(j+j-1) = y(j);
    xe(j+j-1) = x(j);
end

for j = 1:ns
    if mod(j,2) == 1
        Ga_ref = append(Ga_ref,beam_min_list(:,:,length(nFE)));
    elseif mod(j,2) == 0
        Ga_ref = append(Ga_ref,S);
    end
    if mod(j,2) == 1
        if ye(j) == Inf
            PIl = blkdiag(PIl,zeros(pj));
            PIr = blkdiag(PIr,zeros(mj));
            Ga = append(Ga,beam_min_list(:,:,xe(j)));
        else
            PIl = blkdiag(PIl,eye(pj));
            PIr = blkdiag(PIr,eye(mj));
            Ga = append(Ga,ssred(:,:,ye(j),xe(j)));
        end
    elseif mod(j,2) == 0
        PIl = blkdiag(PIl,zeros(psd));
        PIr = blkdiag(PIr,zeros(msd));
        Ga = append(Ga,S);
    end
end



% Creating the closed-loop model that has some of the subsystems reduced
% 'type1'
sscl = rr*feedback(ssdiag,kk,+1)*hh; % same sscl = rr*feedback(ssdiag,kk,-1)*hh;

% Creating the closed-loop model that has some of the subsystems reduced
% 'type2'
Gclr = R*feedback(Ga,K,+1)*H; % same sscl = rr*feedback(ssdiag,kk,-1)*hh;

% % % Testing sscl with
% K = kk; H = hh; R = rr;
% A_test = ssdiag.A + ssdiag.B*K*ssdiag.C;
% B_test = ssdiag.B*H;
% C_test = R*ssdiag.C;
% sscl_test = ss(A_test, B_test, C_test, zeros(4,2));
% sscl_testt = R*(eye(size(ssdiag*K)) - ssdiag*K)^(-1)*ssdiag*H;
% sscl_testtt = R*inv(eye(size(ssdiag*K)) - ssdiag*K)*ssdiag*H;
% bodemag(sscl,'--b',sscl_test,'--g',sscl_testt,'--r',sscl_testtt,'--y');

% Creating the closed-loop reference model 'type1'
cl_ref = rr*feedback(ssdiag_ref,kk,+1)*hh; % same cl_ref = rr*feedback(ssdiag_ref,kk,-1)*hh;

% Creating the closed-loop reference model
Gcl_ref = R*feedback(Ga_ref,K,+1)*H; % same cl_ref = rr*feedback(ssdiag_ref,kk,-1)*hh;

% Closed-loop system error 'type1'
sscl_err = cl_ref - sscl;
sscl_err_norm = norm(sscl_err,'inf');

% Closed-loop system error 'type2'
Gclr_err = Gcl_ref - Gclr;
Gclr_err_norm = norm(Gclr_err,'inf');

% "Gains" of coupled subsystems 'type1'
g1 = norm(pir * kk * feedback(eye(size(pil)), ssdiag * kk, +1),'inf');
g2 = norm(rr * feedback(eye(size(pil)), ssdiag * kk, +1) * pil,'inf');
g3 = norm(feedback(eye(size(pir)), kk * ssdiag, +1) * kk * pil,'inf');
g4 = norm(pir * feedback(eye(size(pir)), kk * ssdiag, +1) * hh,'inf');

% "Gains" of coupled subsystems 'type2'
g1e = norm(PIr * K * feedback(eye(size(PIl)), Ga * K, +1),'inf');
g2e = norm(R * feedback(eye(size(PIl)), Ga * K, +1) * PIl,'inf');
g3e = norm(feedback(eye(size(PIr)), K * Ga, +1) * K * PIl,'inf');
g4e = norm(PIr * feedback(eye(size(PIr)), K * Ga, +1) * H,'inf');

% Indexing gamma corresponding to the x and y vectors
loc = [y',x'];
loc(any(isinf(loc),2),:) = []; % gamma is not calculated for unreduced system so it is discarded
gamma_ind = sub2ind(size(gammau),loc(:,1), loc(:,2));
gamma = max(gammau(gamma_ind));
% gamma_test = 2*max(gammau_test(gamma_ind));

% Checking condition 3.12 (from @REIS2008) 'type1'
if 2*max(g1,g3)*gamma < 1
    fprintf('CLnorm satisfied.\n')
    c_1 = 2 * g2 * (norm(hh, 2) + g1 * norm((ssdiag * hh),'inf'));
    c_2 = 2 * g4 * (norm(rr, 2) + g3 * norm((rr * ssdiag),'inf'));
    % Calculating 3.14 (from @REIS2008)
    gammasscl = min(c_1, c_2) * gamma;
else
    faulty_condition = 2*max(g1,g3)*gamma;
    fprintf('CLnorm NOT satisifed!\n')
    fprintf('faulty_condition = %4.4e \n', faulty_condition)
    c_1 = 2 * g2 * (norm(hh, 2) + g1 * norm((ssdiag * hh),'inf'));
    c_2 = 2 * g4 * (norm(rr, 2) + g3 * norm((rr * ssdiag),'inf'));
    % Calculating 3.14 (from @REIS2008)
    gammasscl = min(c_1, c_2) * gamma;
end

% Checking condition 3.12 (from @REIS2008) 'type2'
if 2*max(g1e,g3e)*gamma < 1
    fprintf('CLnorm satisfied.\n')
    c_1e = 2 * g2e * (norm(H, 2) + g1e * norm((Ga * H),'inf'));
    c_2e = 2 * g4e * (norm(R, 2) + g3e * norm((R * Ga),'inf'));
    % Calculating 3.14 (from @REIS2008)
    gammaGclr = min(c_1e, c_2e) * gamma;
else
    faulty_condition_e = 2*max(g1e,g3e)*gamma;
    fprintf('CLnorm NOT satisifed!\n')
    fprintf('faulty_condition = %4.4e \n', faulty_condition_e)
    c_1e = 2 * g2e * (norm(H, 2) + g1e * norm((Ga * H),'inf'));
    c_2e = 2 * g4e * (norm(R, 2) + g3e * norm((R * Ga),'inf'));
    % Calculating 3.14 (from @REIS2008)
    gammaGclr = min(c_1e, c_2e) * gamma;
end

% Cheking if the condition 3.14 (from @REIS2008) really holds true 'type1'
if sscl_err_norm < gammasscl
    fprintf('Norm of the closed loop subsystem is lower than bound.\n')
    fprintf('sscl_err_norm = %4.4e < %4.4e = gammasscl \n',sscl_err_norm, gammasscl)
else
    fprintf('Norm of the closed loop subsystem is ABOVE the bound!\n')
    fprintf('sscl_err_norm = %4.4e > %4.4e = gammasscl \n',sscl_err_norm, gammasscl)
end

% Cheking if the condition 3.14 (from @REIS2008) really holds true 'type2'
if Gclr_err_norm < gammaGclr
    fprintf('Norm of the closed loop subsystem is lower than bound.\n')
    fprintf('Gclr_err_norm = %4.4e < %4.4e = gammaGclr \n',Gclr_err_norm, gammaGclr)
else
    fprintf('Norm of the closed loop subsystem is ABOVE the bound!\n')
    fprintf('Gclr_err_norm = %4.4e > %4.4e = gammaGclr \n',Gclr_err_norm, gammaGclr)
end

% temp = load('test.txt');
% fileID = fopen('test.txt','a');
% fprintf(fileID, '%1i %2.7e\n',temp(end,1)+1, (gammasscl-sscl_err_norm));
% fclose(fileID);

% for j = 1:nd
%     for jj = 1:msord
%         fid = fopen( sprintf( 'nugap_local%i.txt',j),"A" );
%         [nu, nugap] = gapmetric(beam_min_list(:,:,j),ssred(:,:,jj,j));
%         fprintf(fid, '%1i %2.7e %3.7e\n',jj, nu, nugap);
%         fclose(fid);
%     end
% end

% for j = 1:nd
%     for jj = 1:msord
%         fid = fopen( sprintf( 'nugap_local2ref%i.txt',j),"A" );
%         [nu, nugap] = gapmetric(beam_min_list(:,:,10),ssred(:,:,jj,j));
%         fprintf(fid, '%1i %2.7e %3.7e\n',jj, nu, nugap);
%         fclose(fid);
%     end
% end

% wn = {};
% zeeta = {};
% poul = {};
% for j = 1:nd
%     [wn{j},zeeta{j},poul{j}] = damp(beam_ss_list(:,:,j));
% end
%
% for j = 1:nd
%     wen(j) = max(wn{j});
%     zueta(j) = max(zeeta{j});
%     pouls(j) = max(poul{j});
% end

% % % Plot sigma and lower and upper error for each discretization
% x0=0;
% y0=0;
% width=8.4;
% height=4;
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
% nd_plot = 3; % chose which discretization you want to plot for
% plot(ssredord(1:length(ssredord)),ssrederrnorm(:,nd_plot),'k',ssredord(1:length(ssredord)),gammau(:,nd_plot),'--r',ssredord(1:length(ssredord)),gammal(:,nd_plot),'--g')
% ylabel('$|| \it{G}_{3} - \tilde{G}_{3,o}||$', 'interpreter','latex')
% legend('$|| \it{G}_{3} - \tilde{G}_{3,o}||$','upper bound','lower bound','interpreter','latex')

% % Plot the aboslute error of the second subsystem and the upper bound gamma2
% gamma2TF = tf(gamma2(n_syst),1);
% gamma2TF_4_2 = [gamma2TF, gamma2TF; gamma2TF, gamma2TF; gamma2TF, gamma2TF; gamma2TF, gamma2TF];
% bodemag(ssrederr(:,:,jk,j), 'k', gamma2TF_4_2, '-.k',{1e1,1e6})

% % Plot the aboslute error of the reduced closed-loop subsystem and reference closed-loop subsystem and the accompying upper bounds
% % Minimizing cl_ref
% hsv_tol_cl = 1e-10;
% [~,cl_ref_temp_info] = balred(beam_ss_list(:,:,j),balred_opts);
% beam_hsv = cl_ref_temp_info.HSV; % Hankel singular values
% red_ord_cl = sum(abs(beam_hsv) > hsv_tol_cl);
% cl_ref_temp = balred(cl_ref, red_ord_cl, balred_opts); % minimal order beam state space model
% l = 22;
% [cl_ref_red,cl_bal_info] = balred(cl_ref_temp,l,balred_opts);
% cl_ref_err = cl_ref_red - cl_ref_temp;
% gammacl = 2 * sum(cl_bal_info.HSV(l+1:end));
% gammaclTF = tf(gammacl,1);
% gammaclTF = [gammaclTF, gammaclTF; gammaclTF, gammaclTF; gammaclTF, gammaclTF; gammaclTF, gammaclTF];
% gammassclTF = tf(gammasscl,1);
% gammassclTF = [gammassclTF, gammassclTF; gammassclTF, gammassclTF; gammassclTF, gammassclTF; gammassclTF, gammassclTF];
% x0=0;
% y0=0;
% width=8.4;
% height=4;
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
% bodemag(sscl_err(1,1),'-.r',gammassclTF(1,1),'--r', cl_ref_err(1,1),'k',gammaclTF(1,1),':k',{0.9e+1,1e6}) % Fig 4.2 (b)
% ylabel('$\gamma^{abs} \rm{ and } \gamma$', 'interpreter','latex')
% legend('$\gamma_{cl}^{abs}$', ...
%     '$\gamma_{int}$', ...
%     '$\gamma_{int}^{abs}$', ...
%     '${\gamma}_{cl}$','interpreter','latex')
