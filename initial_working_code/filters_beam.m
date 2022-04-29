if exist('ssred','var') && ~isempty('ssred')
    disp('Beams found. Continuing.')
    return
else
    disp('Creating beams.')
    beams;
end

% Making input filter that covers frequency response per channel
w = logspace(1,5,200);
% sys = frd(beam_minimal_list_sub_ex(:,:,1),w); & krivo, treba oriinalni
% 2x2 sustav
% sys = frd((beam_min_list(1:2,1:2,end) - beam_min_list(1:2,1:2,1)),w); vs nominal low FE
bro = 5;
sys = frd((beam_min_list(1:2,1:2,end) - ssred(1:2,1:2,bro,1)),w); % vs reduced
sysg = abs(sys); % jel mi to treba?
ord = 3;
[no,ni] = iosize(sys);
abc = rss(ord,no,ni); % preallocation
for j = 1:no
    for k = 1:ni
        C.UpperBound = [];
        C.LowerBound = sysg(j,k);
        b = fitmagfrd(sysg(j,k),ord,[],sysg(j,k),C);
        abc(j,k) = b;
    end
end
% bodemag(sysg,':k.',abc,'--g')

% Creating high fidelity reduced order model of the approximated filter to avoid numerical
% complications
hsv_tol_abc = 1e-12;
[~,abc_red_info] = balred(abc,balred_opts);
abc_hsv = abc_red_info.HSV; % Hankel singular values
red_ord_abc = sum(abs(abc_hsv) > hsv_tol_abc);
abc_min = balred(abc, red_ord_abc, balred_opts); % minimal order beam state space model
% bodemag(sysg,':k.',abc_min,'--g')

% M = [zeros(no,ni) abc_min; abc_min beam_min_list(1:2,1:2,1)]; nominal low FE
M = [zeros(no,ni) abc_min; abc_min ssred(1:2,1:2,bro,1)];
M.InputName = {'p(1)' 'p(2)' 'w(1)' 'w(2)'};
M.OutputName = {'q(1)' 'q(2)' 'z(1)' 'z(2)'};
hsv_tol_M = 1e-12;
[~,M_red_info] = balred(M,balred_opts);
M_hsv = M_red_info.HSV; % Hankel singular values
red_ord_M = sum(abs(M_hsv) > hsv_tol_M);
M_min = balred(M, red_ord_M, balred_opts); % minimal order beam state space model
deltoa = ultidyn('delta',[2,2],'Bound',0.1);
deltoa.InputName = 'q';
deltoa.OutputName = 'p';
Mdeltoa = connect(M,deltoa,'w','z');
% bodemag(Mdeltoa,':k.', abc_min, '--r',{1,1e5})


[usys,info_ucover] = ucover(stack(1,beam_min_list(1:2,1:2,end),beam_min_list(1:2,1:2,end-1)),beam_min_list(1:2,1:2,1),[4,4],[],'Additive');
% bodemag(usys,':k.',abc_min,'--r',{1,1e5})

