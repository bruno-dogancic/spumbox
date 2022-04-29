%% Choosing number of systems and their properties (number of finite elements, order of reduction ...)
nb = 5;%10;% % number of beams
ns = nb + nb -1; % number of systems (beams and springs between)

nd = 10;% % number of discretizations
x = [nd, nd, nd, nd, nd];%, nd, nd, 4, 4, 3];%

msord = 20; % maximum order of the reduced subsystems
y = [6, 3, 3, 2, 2];%, 2, 2, 2, 2, 1];%, msord-10, msord, msord, msord-7, msord-13];%
% y = [10, 6, 3, 2, 2];%, 2, 2, 2, 2, 1];%, msord-10, msord, msord, msord-7, msord-13];%

% Vector for choosing orders of W1
z = [4, 3, 2, 2, 2];%, 2, 2, 2, 2, 1];%

% Vector for choosing orders of Wz
zz = [0, 3, 2, 2, 2];%, 2, 2, 2, 1, 1];% % for all systems except first

% Vector for choosing orders of Ww
zw = [0, 3, 2, 2, 2];%, 2, 2, 1, 1, 1];% % for all systems except first

% frequency range of interest
wlog = logspace(1,5,200);

% Dynamical system with forces as inputs and outputs can be:
esins = [1:2]; % forces as (extended system) inputs
% esouts = [1:2]; % only velocities as (extended system) outputs, yields passive systems!
% esouts = [3:4]; % only displacements as (extended system) outputs
esouts = [1:4]; % both velocities and displacements as (extended system) outputs


%% Creating extended subsystem choosing vectors
ye = zeros(1,length(y)+length(y)-1);
xe = zeros(1,length(x)+length(x)-1);
G = struct;
for j = 1:nb
    ye(j+j-1) = y(j);
    xe(j+j-1) = x(j);
    G(j+j-1).W1ord = z(j);
    G(j+j-1).Wzord = zz(j);
    G(j+j-1).Wword = zw(j);
end
tic
%% Series of discretized beams
% [beam_ss_list,~, beam_min_list, ~] = sdb(nd); % series of discretized beams
[beam_ss_list,~, beam_min_lista, ~] = sdb(nd); % series of discretized beams
beam_min_list = beam_min_lista(esouts,esins,:); % extracting subsystem (velocities/displacement)

%% Creating reduced models of subsystems
balred_opts = balredOptions('StateProjection','truncate','ErrorBound','absolute');
for j = 1:nd
    ssredord = 1:msord; % series of orders for model reduction of second system
    [beam_red(:,:,:,j),beam_bal_info(j)] = balred(beam_min_list(:,:,j),ssredord,balred_opts); % series of reduced second subsystems
end
disp('Series of discretized models and BTM reduced order models finished.')
toc
%% Connecting systems
% Dampers' viscous coefficients
c1 = 1e-2; % [Ns/m] [kg/s]
c2 = 1e-1; %

% Springs' stiffness
% k1 = 1e4; % [N/m] [kg/s^2] default
% k2 = 1e5; % default
k1 = 2e2; % [N/m] [kg/s^2]
k2 = 1e1; %

% spring-damper system
sprdam = ss([-c1  0  -k1  0; 0  -c2  0  -k2]);
% for extracting only damping sprdam(:,1:2)



for j = 1:ns
    if mod(j,2) == 1 % odd numbered systems are beams, G1, G3, G5 ...
        if ye(j) == Inf
            G(j).ssr = beam_min_list(:,:,xe(j));
        elseif ye(j) ~= Inf
            G(j).ssr = beam_red(:,:,ye(j),xe(j));
        end
        G(j).ss = beam_min_list(:,:,xe(j));
        [G(j).no,G(j).ni] = iosize(G(j).ss);
        G(j).sdelta = ss(zeros(G(j).ni,G(j).no)); % static
        G(j).sW1 = ss(zeros(G(j).no,G(j).ni)); % static
        G(j).sW2 = ss(zeros(G(j).no,G(j).ni)); % static
        %         G(j).sdelta = ss(eye(G(j).ni,G(j).no)); % static
        %         G(j).sW1 = ss(eye(G(j).no,G(j).ni)); % static
        %         G(j).sW2 = ss(eye(G(j).no,G(j).ni)); % static

    elseif mod(j,2) == 0 % even numbered systems are spring-dampers, G2, G4, G6 ...
        if size(G(j-1).ss.y,1) == 2 % only springs or dampers
            G(j).ss = sprdam(esins,esouts);
            G(j).ssr = sprdam(esins,esouts);
        elseif size(G(j-1).ss.y,1) == 4
            G(j).ss = sprdam;
            G(j).ssr = sprdam;
        else
            error('There is a problem with the number of inputs and outputs of beams and spring-dampers...')
        end
        [G(j).no,G(j).ni] = iosize(G(j).ss);

    end
end

%% Analysis points, summing junctions and branching points are same for all compositions of systems
% Leave it commented before testing
% for j = 1:ns
%     if j == 1
%         % Analysis point on the output of G
%         G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j).no);
%         for jk = 1:length(G(j).apz)
%             G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
%         end
%         % Analysis point on the input of G
%         G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j).ni);
%         for jk = 1:length(G(j).apw)
%             G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
%         end
%         G(j).sum = sumblk(strcat('w',num2str(j),' = w + z',num2str(j+1)),G(j).ni); % summing junction
%         G(j).bp = sumblk(strcat('z = z',num2str(j)),G(j).no); % branching point
%     elseif j >= 2 % || k <= nS - 1
%         if mod(j,2) == 0 % even systems are spring dampers, G2, G4 ...
%             % AP z
%             G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j-1).ni);
%             for jk = 1:length(G(j).apz)
%                 G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
%             end
%             % AP w
%             G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j-1).no);
%             for jk = 1:length(G(j).apw)
%                 G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
%             end
%         elseif mod(j,2) == 1 % odd systems are beams, G3, G5 ...
%             G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j).no);
%             % AP z
%             for jk = 1:length(G(j).apz)
%                 G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
%             end
%             % AP w
%             G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j).ni);
%             for jk = 1:length(G(j).apw)
%                 G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
%             end
%         end
%         if j == 2
%             G(j).sum = sumblk(strcat('w',num2str(j),' = z',num2str(j-1), ' + z',num2str(j+1)),G(j).ni); % summing junction
%         elseif j >= 3 && j <= ns - 1
%             G(j).sum = sumblk(strcat('w',num2str(j),' = z',num2str(j-1), ' + z',num2str(j+1)),G(j).ni); % summing junction
%         elseif j == ns
%             G(j).bp = sumblk(strcat('w',num2str(j),' = z',num2str(j-1)),G(j).ni); % summing junction
%         end
%     end
% end


%% Original system G(end).agg
for j = 1:ns
    if j == 1
        % Analysis point on the output of G
        G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j).no);
        for jk = 1:length(G(j).apz)
            G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
        end
        % Analysis point on the input of G
        G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j).ni);
        for jk = 1:length(G(j).apw)
            G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
        end
        G(j).swdw = G(j).sW2*G(j).sdelta*G(j).sW1;
        G(j).sssb = G(j).ss + G(j).swdw;
        G(j).sys = G(j).apz*G(j).sssb*G(j).apw;
        G(j).sys.InputName = strcat('w',num2str(j));
        G(j).sys.OutputName = strcat('z',num2str(j));

        G(j).sum = sumblk(strcat('w',num2str(j),' = w + z',num2str(j+1)),G(j).ni); % summing junction
        G(j).bp = sumblk(strcat('z = z',num2str(j)),G(j).no); % branching point
        G(j).agg = connect(G(j).sys,G(j).sum,G(j).bp,{'w',strcat('z',num2str(j+1))},{'z',strcat('z',num2str(j))});
    elseif j >= 2 % || k <= nS - 1
        if mod(j,2) == 0 % even systems are spring dampers, G2, G4 ...
            % AP z
            G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j-1).ni);
            for jk = 1:length(G(j).apz)
                G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
            end
            % AP w
            G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j-1).no);
            for jk = 1:length(G(j).apw)
                G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
            end
            G(j).sys = G(j).apz*G(j).ss*G(j).apw;
            G(j).sys.InputName = strcat('w',num2str(j));
            G(j).sys.OutputName = strcat('z',num2str(j));

        elseif mod(j,2) == 1 % odd systems are beams, G3, G5 ...
            G(j).apz = AnalysisPoint(strcat('Z',num2str(j)),G(j).no);
            % AP z
            for jk = 1:length(G(j).apz)
                G(j).apz.Location(jk) = {strcat('z',num2str(j),'(',num2str(jk),')')};
            end
            % AP w
            G(j).apw = AnalysisPoint(strcat('W',num2str(j)),G(j).ni);
            for jk = 1:length(G(j).apw)
                G(j).apw.Location(jk) = {strcat('w',num2str(j),'(',num2str(jk),')')};
            end
            G(j).swdw = G(j).sW2*G(j).sdelta*G(j).sW1;
            G(j).sssb = G(j).ss + G(j).swdw;
            G(j).sys = G(j).apz*G(j).sssb*G(j).apw;
            G(j).sys.InputName = strcat('w',num2str(j));
            G(j).sys.OutputName = strcat('z',num2str(j));
        end
        if j == 2
            G(j).sum = sumblk(strcat('w',num2str(j),' = z',num2str(j-1), ' + z',num2str(j+1)),G(j).ni); % summing junction
            G(j).agg = connect(G(j-1).agg,G(j).sys,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j >= 3 && j <= ns - 1
            G(j).sum = sumblk(strcat('w',num2str(j),' = z',num2str(j-1), ' + z',num2str(j+1)),G(j).ni); % summing junction
            G(j).agg = connect(G(j-1).agg,G(j).sys,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});%, strcat('z',num2str(k-1)) ); % we add analysis point from previous system
        elseif j == ns
            G(j).bp = sumblk(strcat('w',num2str(j),' = z',num2str(j-1)),G(j).ni); % branching point
            G(j).agg = connect(G(j-1).agg,G(j).sys,G(j).bp,'w', ...
                'z');%,{'z2','z3','z4','z5','z6','z7'});%, {strcat('z',num2str(k-1)), strcat('z',num2str(k))} ); % we add analysis point from previous and last system
        end
    end
end

%% Generating weighting filters W1 and W2 and unstructured uncertainty |delta| <= 1

for j = 1:ns
    if mod(j,2) == 1 % only for odd numbered systems (beams) G1, G3, G5 ...
        G(j).delta = ultidyn(strcat('delta',num2str(j)),[G(j).ni,G(j).no]);
        G(j).rel = G(j).ss - G(j).ssr;

        % W1
        G(j).W1frd = frd(G(j).rel,wlog);
        G(j).W1frda = abs(G(j).W1frd);
        G(j).W1 = rss(G(j).W1ord,G(j).no,G(j).ni); % preallocation of filter W1 approximating connected enviorement
        for jj = 1:G(j).no
            for kk = 1:G(j).ni
                C1.UpperBound = [];
                C1.LowerBound = G(j).W1frda(jj,kk);
                b1 = fitmagfrd(G(j).W1frda(jj,kk),G(j).W1ord,[],G(j).W1frda(jj,kk),C1);
                G(j).W1(jj,kk) = b1;
            end
        end
        %         % W2, setup for ucover
        %         G(j).W2frd = frd(G(j).rel,wlog);
        %         G(j).W2frda = abs(G(j).W1frd);
        %         G(j).W2ord = 3; % order of filter W11 that approximates additive model G1 - G1t
        %         G(j).W2 = rss(G(j).W1ord,G(j).no,G(j).ni); % preallocation of filter W1 approximating connected enviorement
        %         for jj = 1:G(j).no
        %             for kk = 1:G(j).ni
        %                 G(jj).C2.UpperBound = [];
        %                 G(jj).C2.LowerBound = G(j).W2frda(jj,kk);
        %                 G(jj).b2 = fitmagfrd(G(j).W2frda(jj,kk),G(j).W2ord,[],G(j).W2frda(jj,kk),G(jj).C2);
        %                 G(j).W2(jj,kk) = G(jj).b2;
        %             end
        %         end
        % W2 can be used to further decrease uncertainty, or together with W1 to
        % get two smaller order filters - for now we use only identity W2
        G(j).W2 = ss(eye(G(j).no,G(j).ni));
    end
end

%% Creating reduced order aggregate model, with previously obtained filters
% We don't need analysis points anymore and use same summing
% junctions and branching points as before
for j = 1:ns
    if j == 1
        G(j).wdw = G(j).W2*G(j).delta*G(j).W1;
        G(j).sysr = G(j).ssr + G(j).wdw;
        G(j).sysr.InputName = strcat('w',num2str(j));
        G(j).sysr.OutputName = strcat('z',num2str(j));
        G(j).aggr = connect(G(j).sysr,G(j).sum,G(j).bp,{'w',strcat('z',num2str(j+1))},{'z',strcat('z',num2str(j))});
    elseif j >= 2 % && j <= ns - 1
        if mod(j,2) == 0 % even systems are spring dampers, G2, G4 ...
            G(j).sysr = G(j).ssr;
            G(j).sysr.InputName = strcat('w',num2str(j));
            G(j).sysr.OutputName = strcat('z',num2str(j));
        elseif mod(j,2) == 1 % odd systems are beams, G3, G5 ...
            G(j).wdw = G(j).W2*G(j).delta*G(j).W1;
            G(j).sysr = G(j).ssr + G(j).wdw;
            G(j).sysr.InputName = strcat('w',num2str(j));
            G(j).sysr.OutputName = strcat('z',num2str(j));
        end
        if j == 2
            G(j).aggr = connect(G(j-1).aggr,G(j).sysr,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j >= 3 && j <= ns - 1
            G(j).aggr = connect(G(j-1).aggr,G(j).sysr,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j == ns
            G(j).aggr = connect(G(j-1).aggr,G(j).sysr,G(j).bp,'w', ...
                'z');
        end
    end
end

%% Uncertainty (weighting filters) refinement
% Getting input-output transfer function replacing the surrounding
for j = 1:ns
    if mod(j,2) == 1
        if j >= 3 && j < ns-1
            G(j).(strcat('Wz',num2str(j))) = ss(getIOTransfer(G(end).agg,strcat('z',num2str(j)),'z',strcat('z',num2str(j))));
            G(j).(strcat('Wz',num2str(j))).InputName = strcat('z',num2str(j));
            G(j).(strcat('Wz',num2str(j))).OutputName = 'z';

            G(j).(strcat('Wz',num2str(j-1))) = ss(getIOTransfer(G(end).agg,'w',strcat('z',num2str(j-1)),strcat('z',num2str(j))));
            G(j).(strcat('Wz',num2str(j-1))).InputName = 'w';
            G(j).(strcat('Wz',num2str(j-1))).OutputName = strcat('z',num2str(j-1));

            %         G(j).(strcat('Wz',num2str(j+1))) = ss(getIOTransfer(G(end).agg,strcat('z',num2str(j)),strcat('z',num2str(j+1)),{strcat('z',num2str(j)),strcat('z',num2str(j+1))}));
            G(j).(strcat('Wz',num2str(j+1))) = ss(getIOTransfer(G(end).agg,strcat('z',num2str(j)),strcat('z',num2str(j+1)),strcat('z',num2str(j))));
            G(j).(strcat('Wz',num2str(j+1))).InputName = strcat('z',num2str(j));
            G(j).(strcat('Wz',num2str(j+1))).OutputName = strcat('z',num2str(j+1));
        elseif j == ns
            G(j).(strcat('Wz',num2str(j))) = ss(getIOTransfer(G(end).agg,strcat('z',num2str(j)),'z',{strcat('z',num2str(j)),strcat('w',num2str(j))}));
            G(j).(strcat('Wz',num2str(j))).InputName = strcat('z',num2str(j));
            G(j).(strcat('Wz',num2str(j))).OutputName = 'z';

            G(j).(strcat('Wz',num2str(j-1))) = ss(getIOTransfer(G(end).agg,'w',strcat('z',num2str(j-1)),{strcat('z',num2str(j)),strcat('w',num2str(j))}));
            G(j).(strcat('Wz',num2str(j-1))).InputName = 'w';
            G(j).(strcat('Wz',num2str(j-1))).OutputName = strcat('z',num2str(j-1));
        end
    end
end

% Creating input and output scaling weights - low order filter from
% previously obtained input output transfer functions

% Wz are filters of reduced surrounding of the system at each (j >=3 to
% last) system output, hence z in name


% Ww are filters of reduced surrounding of the system at each (j >=3 to
% last) system input, hence w in name


for j = 1:ns
    if mod(j,2) == 1
        if j >= 3 && j < ns-1
            G(j).Wzfrd = frd(G(j).(strcat('Wz',num2str(j))),wlog); % zj --> z, Wzj
            G(j).Wzfrda = abs(G(j).Wzfrd);
            G(j).Wz = rss(G(j).Wzord,G(j).no,G(j).ni);
            for jj = 1:G(j).no
                for kk = 1:G(j).ni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).Wzfrda(jj,kk);
                    b1 = fitmagfrd(G(j).Wzfrda(jj,kk),G(j).Wzord,[],G(j).Wzfrda(jj,kk),C1);
                    G(j).Wz(jj,kk) = b1;
                end
            end
            G(j).Wz.InputName = strcat('z',num2str(j));
            G(j).Wz.OutputName = 'z';

            % At this moment I'm not sure if I should include dynamics of the
            % systems below j-th system for uncertainty (filters) refinement or should I
            % only take the transfer function form external input w -->
            % z(j-1)????
            % So basically if I want to do the same thing I did in
            % beam_connections.m I have to use only upper part and put
            % Below is such case:
            %%%%%
            % G(j).Wwfrd = frd(G(j).(strcat('Wz',num2str(j-1))),wlog); % w --> z(j-1), Wz(j-1)
            %%%%%
            % I will proceed for now with using both TF and summing them at G(j).sum (or parallel systems):
            %%%%%
            % G(j).Wwfrd = frd( parallel(G(j).(strcat('Wz',num2str(j-1))), ...
            %    G(j).(strcat('Wz',num2str(j+1)))), wlog); %sum of w --> z(j-1), Wz(j-1) and z(j) --> z(j+1), Wz(j+1)
            %%%%%
            % I have come to conclusion that it is best to use low order of
            % both and connect them, such that the lower acts as feedback,
            % and upper as an external input to the system. For that I will
            % rename them accordingly, upper will be w1 and lower will be
            % w2 (all will use same order Wword)
            G(j).Ww1frd = frd(G(j).(strcat('Wz',num2str(j-1))),wlog); % w --> z(j-1), Wz(j-1)
            G(j).Ww1frda = abs(G(j).Ww1frd);
            [tempno,tempni] = iosize(G(j).(strcat('Wz',num2str(j-1))));
            G(j).Ww1 = rss(G(j).Wword,tempno,tempni);
            for jj = 1:tempno
                for kk = 1:tempni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).Ww1frda(jj,kk);
                    b1 = fitmagfrd(G(j).Ww1frda(jj,kk),G(j).Wword,[],G(j).Ww1frda(jj,kk),C1);
                    G(j).Ww1(jj,kk) = b1;
                end
            end
            G(j).Ww1.InputName = 'w';
            G(j).Ww1.OutputName = strcat('z',num2str(j-1));

            G(j).Ww2frd = frd(G(j).(strcat('Wz',num2str(j+1))),wlog); % zj --> z(j+1), Wz(j+1)
            G(j).Ww2frda = abs(G(j).Ww2frd);
            [tempno,tempni] = iosize(G(j).(strcat('Wz',num2str(j+1))));
            G(j).Ww2 = rss(G(j).Wword,tempno,tempni);
            for jj = 1:tempno
                for kk = 1:tempni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).Ww2frda(jj,kk);
                    b1 = fitmagfrd(G(j).Ww2frda(jj,kk),G(j).Wword,[],G(j).Ww2frda(jj,kk),C1);
                    G(j).Ww2(jj,kk) = b1;
                end
            end
            G(j).Ww2.InputName = strcat('z',num2str(j));
            G(j).Ww2.OutputName = strcat('z',num2str(j+1));

        elseif j == ns

            G(j).Wzfrd = frd(G(j).(strcat('Wz',num2str(j))),wlog); % zj --> z, Wzj
            G(j).Wzfrda = abs(G(j).Wzfrd);
            [tempno,tempni] = iosize(G(j).(strcat('Wz',num2str(j))));
            G(j).Wz = rss(G(j).Wzord,tempno,tempni);
            for jj = 1:tempno
                for kk = 1:tempni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).Wzfrda(jj,kk);
                    b1 = fitmagfrd(G(j).Wzfrda(jj,kk),G(j).Wzord,[],G(j).Wzfrda(jj,kk),C1);
                    G(j).Wz(jj,kk) = b1;
                end
            end

            % For this case it is easy, since there is nothing below the last
            % system :)
            G(j).Wwfrd = frd(G(j).(strcat('Wz',num2str(j-1))),wlog); % w --> z(j-1), Wz(j-1)
            G(j).Wwfrda = abs(G(j).Wwfrd);
            [tempno,tempni] = iosize(G(j).(strcat('Wz',num2str(j-1))));
            G(j).Ww = rss(G(j).Wword,tempno,tempni);
            for jj = 1:tempno
                for kk = 1:tempni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).Wwfrda(jj,kk);
                    b1 = fitmagfrd(G(j).Wwfrda(jj,kk),G(j).Wword,[],G(j).Wwfrda(jj,kk),C1);
                    G(j).Ww(jj,kk) = b1;
                end
            end
        end
    end
end

% Refining weighting filters at uncertainties
% Names will contain "e" as a reminder of using "envoirment" or as
% "extended" systems with surrounding
% I will keep the same order G(j).W1ord of filters for best comparison
for j = 1:ns
    if mod(j,2) == 1 % only for odd numbered systems (beams) G1, G3, G5 ...
        if j == 1 % for the first system, for now, everything remains the same as before
            G(j).rele = G(j).ss - G(j).ssr;

            % W1
            G(j).W1efrd = frd(G(j).rele,wlog);
            G(j).W1efrda = abs(G(j).W1efrd);
            G(j).W1e = rss(G(j).W1ord,G(j).no,G(j).ni); % preallocation of filter W1 approximating connected enviorement
            for jj = 1:G(j).no
                for kk = 1:G(j).ni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).W1efrda(jj,kk);
                    b1 = fitmagfrd(G(j).W1efrda(jj,kk),G(j).W1ord,[],G(j).W1efrda(jj,kk),C1);
                    G(j).W1e(jj,kk) = b1;
                end
            end
            G(j).W2e = ss(eye(G(j).no,G(j).ni));
        elseif j >= 3 && j < ns-1
            G(j).ssr.InputName = strcat('w',num2str(j));
            G(j).ssr.OutputName = strcat('z',num2str(j));
            G(j).ss.InputName = strcat('w',num2str(j));
            G(j).ss.OutputName = strcat('z',num2str(j));
            G(j).sume = sumblk(strcat('w',num2str(j),'=z',num2str(j-1),'+z',num2str(j+1)),G(j).ni);
            G(j).ssre = connect(G(j).Wz,G(j).ssr,G(j).Ww1,G(j).Ww2,G(j).sume,'w','z');
            G(j).sse = connect(G(j).Wz,G(j).ss,G(j).Ww1,G(j).Ww2,G(j).sume,'w','z');
            G(j).rele = G(j).sse - G(j).ssre;

            % W1
            G(j).W1efrd = frd(G(j).rele,wlog);
            G(j).W1efrda = abs(G(j).W1efrd);
            G(j).W1e = rss(G(j).W1ord,G(j).no,G(j).ni); % preallocation of filter W1 approximating connected enviorement
            for jj = 1:G(j).no
                for kk = 1:G(j).ni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).W1efrda(jj,kk);
                    b1 = fitmagfrd(G(j).W1efrda(jj,kk),G(j).W1ord,[],G(j).W1efrda(jj,kk),C1);
                    G(j).W1e(jj,kk) = b1;
                end
            end
            G(j).W2e = ss(eye(G(j).no,G(j).ni));
        elseif j == ns
            G(j).ssre = G(j).Wz*G(j).ssr*G(j).Ww;
            G(j).ssre.InputName = 'w';
            G(j).ssre.OutputName = 'z';
            G(j).sse = G(j).Wz*G(j).ss*G(j).Ww;
            G(j).sse.InputName = 'w';
            G(j).sse.OutputName ='z';
            G(j).rele = G(j).sse - G(j).ssre;

            % W1
            G(j).W1efrd = frd(G(j).rele,wlog);
            G(j).W1efrda = abs(G(j).W1efrd);
            G(j).W1e = rss(G(j).W1ord,G(j).no,G(j).ni); % preallocation of filter W1 approximating connected enviorement
            for jj = 1:G(j).no
                for kk = 1:G(j).ni
                    C1.UpperBound = [];
                    C1.LowerBound = G(j).W1efrda(jj,kk);
                    b1 = fitmagfrd(G(j).W1efrda(jj,kk),G(j).W1ord,[],G(j).W1efrda(jj,kk),C1);
                    G(j).W1e(jj,kk) = b1;
                end
            end
            G(j).W2e = ss(eye(G(j).no,G(j).ni));
        end
    end
end

%% Creating aggregate model with refined filters
for j = 1:ns
    if j == 1
        G(j).wdwe = G(j).W2*G(j).delta*G(j).W1;
        G(j).sysre = G(j).ssr + G(j).wdwe;
        G(j).sysre.InputName = strcat('w',num2str(j));
        G(j).sysre.OutputName = strcat('z',num2str(j));
        G(j).aggre = connect(G(j).sysre,G(j).sum,G(j).bp,{'w',strcat('z',num2str(j+1))},{'z',strcat('z',num2str(j))});
    elseif j >= 2 % && j <= ns - 1
        if mod(j,2) == 0 % even systems are spring dampers, G2, G4 ...
            G(j).sysre = G(j).ssr;
            G(j).sysre.InputName = strcat('w',num2str(j));
            G(j).sysre.OutputName = strcat('z',num2str(j));
        elseif mod(j,2) == 1 % odd systems are beams, G3, G5 ...
            G(j).wdwe = G(j).W2e*G(j).delta*G(j).W1e;
            G(j).sysre = G(j).ssr + G(j).wdwe;
            G(j).sysre.InputName = strcat('w',num2str(j));
            G(j).sysre.OutputName = strcat('z',num2str(j));
        end
        if j == 2
            G(j).aggre = connect(G(j-1).aggre,G(j).sysre,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j >= 3 && j <= ns - 1
            G(j).aggre = connect(G(j-1).aggre,G(j).sysre,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j == ns
            G(j).aggre = connect(G(j-1).aggre,G(j).sysre,G(j).bp,'w', ...
                'z');
        end
    end
end

%% Creating static filters, covering the peak of the refined uncertainty filter
for j = 1:ns
    if mod(j,2) == 1 && j >= 3
        for jj = 1:G(j).no
            for kk = 1:G(j).ni
                W1static(jj,kk) = norm(G(j).W1e(jj,kk),Inf);
            end
        end
        G(j).W1s = ss(W1static);
        G(j).W2s = ss(eye(G(j).no,G(j).ni));
     end
end

%% Creating aggregate model with static filters
for j = 1:ns
    if j == 1
        G(j).wdws = G(j).W2*G(j).delta*G(j).W1;
        G(j).sysrs = G(j).ssr + G(j).wdws;
        G(j).sysrs.InputName = strcat('w',num2str(j));
        G(j).sysrs.OutputName = strcat('z',num2str(j));
        G(j).aggrs = connect(G(j).sysrs,G(j).sum,G(j).bp,{'w',strcat('z',num2str(j+1))},{'z',strcat('z',num2str(j))});
    elseif j >= 2 % && j <= ns - 1
        if mod(j,2) == 0 % even systems are spring dampers, G2, G4 ...
            G(j).sysrs = G(j).ssr;
            G(j).sysrs.InputName = strcat('w',num2str(j));
            G(j).sysrs.OutputName = strcat('z',num2str(j));
        elseif mod(j,2) == 1 % odd systems are beams, G3, G5 ...
            G(j).wdws = G(j).W2s*G(j).delta*G(j).W1s;
            G(j).sysrs = G(j).ssr + G(j).wdws;
            G(j).sysrs.InputName = strcat('w',num2str(j));
            G(j).sysrs.OutputName = strcat('z',num2str(j));
        end
        if j == 2
            G(j).aggrs = connect(G(j-1).aggrs,G(j).sysrs,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j >= 3 && j <= ns - 1
            G(j).aggrs = connect(G(j-1).aggrs,G(j).sysrs,G(j).sum,{'w',strcat('z',num2str(j+1))}, ...
                {'z',strcat('z',num2str(j))});
        elseif j == ns
            G(j).aggrs = connect(G(j-1).aggrs,G(j).sysrs,G(j).bp,'w', ...
                'z');
        end
    end
end

%% Variables to remove
% % Fields to remove
% flds = {'C1','b1','C1z','b1z','C1w1','b1w1','C1w2','b1w2'};
% G = rmfield(G,flds);
% % Temporary variables to remove
% clear C1 b1 tempno tempni


%% Testing zone

% sprdam = ss([-c1  0  -k1  0;
%              0  -c2  0  -k2]);
% % For testing
% abc = ss(getIOTransfer(G(end).agg,'w','z2','z3'));
% abc_test = (sprdam*G1)*feedback(eye(2,2),sprdam*G1,+1);
% % bodemag(abc,'-k',abc_test,'md')
% % norm(abc,Inf),norm(abc_test,Inf)
%
% def = ss(getIOTransfer(G(end).agg,'w2','z','z3'));
% def_test = (G1*sprdam)*feedback(eye(4,4),G1*sprdam,+1);
% % bodemag(def,'-k',def_test,'md')
% % norm(def,Inf),norm(def_test,Inf)
%
% %
% abc1 = ss(getIOTransfer(G(end).agg,'w3','z3',{'z2','z3'}));
% % bodemag(abc1,'-k',G1,'md')
% % norm(G1,Inf), norm(abc1,Inf)

% checking memory size of G
% memsizeG = whos('G');
% memsizeGmb = memsizeG.bytes*10^-6;

%% Plotting area
% for jf = 3:2:ns
% figure()
% bodemag(G(jf).W1e,':k',G(jf).W1,':r')
% end

%% Generating open loop systems with lftdata
M = lftdata(G(end).agg);
Mr = lftdata(G(end).aggr);
Mre = lftdata(G(end).aggre);
Mrs = lftdata(G(end).aggrs);

%% IQC analysis
% Making indeksing vectors depending on the number of inputs and outputs of
% delta
if G(1).no == 2
    du = zeros(2,nb);
    du(1,:) = 1:G(1).no:(nb*G(1).no-1);
    du(2,:) = G(1).no:G(1).no:nb*G(1).no;
elseif G(1).no == 4
    du = zeros(2,nb);
    du(1,:) = 1:G(1).no:(nb*G(1).no-1);
    du(2,:) = G(1).no:G(1).no:nb*G(1).no;
end
dy = zeros(2,nb);
dy(1,:) = 1:G(1).ni:(nb*G(1).ni-1);
dy(2,:) = G(1).ni:G(1).ni:nb*G(1).ni;

due = zeros(2,length(du)+length(du)-1);
dye = zeros(2,length(dy)+length(dy)-1);
for kk = 1:2
    for jj = 1:nb
        due(kk,jj+jj-1) = du(kk,jj);
        dye(kk,jj+jj-1) = dy(kk,jj);
    end
end
% Creating IQC delta's

clear de pe probr probre probrs

for j=1:ns
    if mod(j,2) == 1
        G(j).de = iqcdelta(strcat('de',num2str(j)),'InputChannel', ...
            [dye(1,j):dye(2,j)],'OutputChannel',[due(1,j):due(2,j)], ...
            'StaticDynamic','D','NormBounds',1,'Structure','FB');
    end
end

de = blkdiag('de',G(1).de,G(3).de);
if nb > 2
    for j = 5:ns
        if mod(j,2) == 1
            de = blkdiag('de',de,G(j).de);
        end
    end
end

de = iqcassign(de,'ultid','Length',3,'PoleLocation',-1, 'BasisFunctionType',1);

pe = iqcdelta('pe','ChannelClass','P','InputChannel',[dye(2,end)+1:dye(2,end)+G(1).ni], ...
    'OutputChannel',[due(2,end)+1:due(2,end)+G(1).no],'PerfMetric','L2');
% tic
% probr = iqcanalysis(Mr,{de,pe},'SolChk ','on','eps',1e-8);
%         probr = iqcanalysis(Mr,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','sdpt3');
% toc

% tic
% probre = iqcanalysis(Mre,{de,pe},'SolChk','on','eps',1e-8);
% % probre = iqcanalysis(Mre,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
% toc

% tic
% probrs = iqcanalysis(Mrs,{de,pe},'SolChk','on','eps',1e-8);
% % probrs = iqcanalysis(Mrs,{de,pe},'SolChk','on','eps',1e-8,'Parser','Yalmip','Solver','mosek');
% toc

% tic
% [~,NUGAP] = gapmetric(G(end).agg,G(end).aggrs)
% toc

% %% Plotting area
% for j = 1:ns
%     if j >= 3
%         if mod(j,2) == 1
%             figure()
%             bodemag(G(j).W1,'-k',G(j).W1e,'--b',G(j).W1s,'-.r')
%         end
%     end
% end

