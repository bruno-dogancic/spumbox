
nd = 10;% % number of d5iscretizations
msord = 23; % maximum order of the reduced subsystems (if it cannot complete check hsv_tol below)

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
    mat_b1 = zeros(2 * nFE(j), 2); % without uncertain channel
    mat_b1(2 * nke3i, 1) = 1;
    mat_b1(4 * nke3i, 2) = 1;


    % Output vector is [v1(t); v2(t); x1(t); x2(t)] - velocities and displacements
    mat_c1 = zeros([2, 2 * nFE(j)]); % without uncertain channel
    mat_c1(1, 2 * nke3i) = 1;
    mat_c1(2, 4 * nke3i) = 1;

    % Creating descriptor state-space dynamical system
    [mat_a_dss, mat_b_dss, mat_c_dss, mat_d_dss, mat_e_dss] = mehdss(mat_m, mat_p, mat_k, mat_b1, mat_c1);
    beam_dss_list(:,:,j) = dss(mat_a_dss,mat_b_dss,mat_c_dss,mat_d_dss,mat_e_dss);

    
    % Creating state-space dynamical system
    [mat_a, mat_b, mat_c, mat_d] = mehss(mat_m, mat_p, mat_k, mat_b1, mat_c1);
    beam_ss_list(:,:,j) = ss(mat_a, mat_b, mat_c, mat_d);
%     beam_ss_list(:,:,j) = ss(beam_dss_list(:,:,j),'explicit'); %  alternative way to calculate state-space
    
    % Creating extended state space beam model with uncertain chanels
%     mat_b_dss_ex = [mat_b_dss, zeros(size(mat_b_dss))]; % adding two inputs
%     mat_b_dss_ex(2 * nke3i - nke3i, 3) = 1; % uncertain force left to the original force at l/3
%     mat_b_dss_ex(4 * nke3i + 2*nke3i, 4) = 1; % uncertain force right to the original force at 2l/3
%     mat_c_dss_ex = [mat_c_dss; zeros(2,size(mat_c_dss,2))]; % adding two outputs
%     mat_c_dss_ex(5, 2 * nke3i - nke3i) = 1; % uncertain displacement left to the original displacement at l/3
%     mat_c_dss_ex(6, 4 * nke3i + 2*nke3i) = 1; % uncertain displacement right to the original displacement at 2l/3
%     mat_d_dss_ex = zeros(size(mat_c_dss_ex,1),size(mat_b_dss_ex,2));
%     beam_dss_list_ex(:,:,j) = dss(mat_a_dss,mat_b_dss_ex,mat_c_dss_ex,mat_d_dss_ex,mat_e_dss);
%% This code bellow (till the "end") is sort of dumb, my extened robust stability and performance model is obtained by interconnecion and obtained from LFT data...
%ne = 1;
%labeltype = 1;
% [beam_dss_list_ex(:,:,j),np,nq,nw,nz] = pqwz(beam_dss_list(:,:,j),ne,labeltype);
% beam_ss_list_ex(:,:,j) = ss(beam_dss_list_ex(:,:,j),'explicit');
% Extracted first two outputs from the overall system - what is the output
% vector?!?!?
% [beam_ss_list_sub_ex(:,:,j),np,nq,nw,nz] = pqwz(beam_ss_list([1,2],:,j),ne,labeltype);
end



%% Minimizing models
% Alocation for each variable that chages size in for loops in this section
% beam_ss_min_list(:,:,length(nke)) = greda_ss_list(:,:,length(nke));

% Creating treshold for truncating small hankel singular values in order to
% get as close to minimal realization as posible, but to leave some space
% for testing the algorithm
hsv_tol = 1e-12; % treshold for truncation
balred_opts = balredOptions('StateProjection','truncate','ErrorBound','absolute');

% Creating series of initially ("high fidelity") reduced discretized models
% which will be used as minimal reference model afterwards
for j = 1:length(nFE)
    [~,beam_red_info] = balred(beam_ss_list(:,:,j),balred_opts);
    beam_hsv = beam_red_info.HSV; % Hankel singular values
    red_ord = sum(abs(beam_hsv) > hsv_tol);
    beam_min_list(:,:,j) = balred(beam_ss_list(:,:,j), red_ord, balred_opts); % minimal order beam state space model
    beam_minimal_list(:,:,j) = minreal(beam_min_list(:,:,j)); % checking if any more states can be discarded
end

% % Creating series of initially ("high fidelity") reduced discretized models
% % which will be used as minimal reference model afterwards
% for j = 1:length(nFE)
%     [~,beam_red_info_ex] = balred(beam_ss_list_ex(:,:,j),balred_opts);
%     beam_hsv_ex = beam_red_info_ex.HSV; % Hankel singular values
%     red_ord_ex = sum(abs(beam_hsv) > hsv_tol);
%     beam_min_list_ex(:,:,j) = balred(beam_ss_list_ex(:,:,j), red_ord_ex, balred_opts); % minimal order beam state space model
%     beam_minimal_list_ex(:,:,j) = minreal(beam_min_list_ex(:,:,j)); % checking if any more states can be discarded
% end

% Creating series of initially ("high fidelity") reduced extended subsystem
% (only first two outputs are considered and extended)
% discretized models which will be used as minimal reference model afterwards
% for j = 1:length(nFE)
%     [~,beam_red_info_sub_ex] = balred(beam_ss_list_sub_ex([1,2],:,j),balred_opts);
%     beam_hsv_sub_ex = beam_red_info_sub_ex.HSV; % Hankel singular values
%     red_ord_ex = sum(abs(beam_hsv_sub_ex) > hsv_tol);
%     beam_min_list_sub_ex(:,:,j) = balred(beam_ss_list_sub_ex(:,:,j), red_ord_ex, balred_opts); % minimal order beam state space model
%     beam_minimal_list_sub_ex(:,:,j) = minreal(beam_min_list_sub_ex(:,:,j)); % checking if any more states can be discarded
% end

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
        gammal(jk,j) = bal_info(j).HSV(jk+1); % for plotting only, no sense since first order reduced system doesn't have lower bound...
    end
end

% %% Creating reduced models of subsystems
% for j = 1:length(nFE)
%     maxrst = msord + 1; % max order of the system
%     ssredord = 1:maxrst-1; % series of orders for model reduction of second system
%     [ssred(:,:,:,j),bal_info(j)] = balred(beam_min_list_sub_ex(:,:,j),ssredord,balred_opts); % series of reduced second subsystems
% 
%     % H_inf norm of subsystem error, upper and lower bounds
%     for jk = 1:length(ssredord)
%         ssrederr(:,:,jk,j) = beam_min_list_sub_ex(:,:,j) - ssred(:,:,jk,j);
%         ssrederrnorm(jk,j) = norm(ssrederr(:,:,jk,j),"inf");
%         gammau(jk,j) = 2 * sum(bal_info(j).HSV(jk+1:end));
%         gammal(jk,j) = bal_info(j).HSV(jk+1); % for plotting only, no sense since first order reduced system doesn't have lower bound...
%     end
% end