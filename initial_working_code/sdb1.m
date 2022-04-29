function [beam_ss_list,beam_dss_list, beam_min_list, beam_minimal_list] = sdb(nd)

%% Series of discretized beams
% For provided input of discretizations, system returns a state-space and
% descriptor state space models that has from 9 with up to 9*nd finite
% elements per model.
% Euler beams are modeled such that inputs to beams are forces at 1/3 and
% 2/3 of the beam's length and outputs are velocities and displacements at
% 1/3 and 2/3 of the beam's length.
%
% Usage:
% [beam_ss_list,beam_dss_list, beam_min_list, beam_minimal_list] = sdb(nd)

e_beam = 210e9; % [GPa] [GN/m^2] [Gkg/s^2/m]
l_beam = 2; % [m]
d_beam = 0.018; % [m]
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
    [~, eigfreq] = eig(mat_k, mat_m,'vector');
    omega = sort(sqrt(eigfreq));

    % damping matrix
    zeta = 0.08; % modal damping coefficient
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

