load('highOrderModel.mat','G')
% rng(1,'twister'); % For reproducibility
% G = rss(15,1,1);
% G = ss(tf([1 0.1 7.5],[1 0.12 9 0 0]));

w = logspace(1e-1,5,200); % frequency range of interest
W = abs(frd(G,w)); % frequency response data
ord = [2,6,8]; % orders of filter W1
for j = 1:length(ord)
  C1.UpperBound = [];
  C1.LowerBound = W;
  W1(j).fit = fitmagfrd(W,ord(j),[],[],C1);
end
w_plot = {1e0,1e2};
bodemag(G,'k',W1(1).fit,'--r',W1(2).fit,'b-.',W1(3).fit,':m',w_plot)
legend('Original system (15th order)','2nd order cover fit','6th order cover fit',...
       '8th order cover fit')