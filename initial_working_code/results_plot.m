%% bc1.m
% x0=10;
% y0=10;
% width=8;
% height=8;
h = figure();
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
bodemag(M('z','w'),'k',Mr('z','w'),'-.r',Mre('z','w'),'--b',Mrs('z','w'),'g',{1e1,1e5})
hax = get(gcf,'Children');
for ii = 1:numel(hax)
if strcmp(hax(ii).Type,'axes')
hc = hax(ii).Children;
for jj = 1:numel(hc)
if strcmp(hc(jj).Type,'hggroup')
if jj == 1
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'o';
hc(jj).Children.MarkerIndices = 1:10:numel(hc(jj).Children.MarkerIndices);
end
if jj == 2
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'x';
hc(jj).Children.MarkerIndices = 2:9:numel(hc(jj).Children.MarkerIndices);
end
if jj == 3
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'd';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
if jj == 4
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'none';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
end
end
end
end

%% bc2.m
% x0=10;
% y0=10;
% width=8;
% height=8;
h = figure();
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
bodemag(M('z','w'),'k',Mr('z','w'),'-.r',Mre('z','w'),'--b',Mrs('z','w'),'g',{1e1,1e5})
hax = get(gcf,'Children');
for ii = 1:numel(hax)
if strcmp(hax(ii).Type,'axes')
hc = hax(ii).Children;
for jj = 1:numel(hc)
if strcmp(hc(jj).Type,'hggroup')
if jj == 1
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'o';
hc(jj).Children.MarkerIndices = 1:10:numel(hc(jj).Children.MarkerIndices);
end
if jj == 2
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'x';
hc(jj).Children.MarkerIndices = 2:9:numel(hc(jj).Children.MarkerIndices);
end
if jj == 3
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'd';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
if jj == 4
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'none';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
end
end
end
end


%% bc3.m
% x0=10;
% y0=10;
% width=8;
% height=8;
h = figure();
% set(gcf,'units','centimeters','position',[x0,y0,width,height])
bodemag(M('z','w'),'k',Mr('z','w'),'-.r',Mre('z','w'),'--b',Mrs('z','w'),'g',{1e1,1e5})
hax = get(gcf,'Children');
for ii = 1:numel(hax)
if strcmp(hax(ii).Type,'axes')
hc = hax(ii).Children;
for jj = 1:numel(hc)
if strcmp(hc(jj).Type,'hggroup')
if jj == 1
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'o';
hc(jj).Children.MarkerIndices = 1:10:numel(hc(jj).Children.MarkerIndices);
end
if jj == 2
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'x';
hc(jj).Children.MarkerIndices = 2:9:numel(hc(jj).Children.MarkerIndices);
end
if jj == 3
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'd';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
if jj == 4
% all Line poperties accessible here
% set Marker and space ever 10th data point
hc(jj).Children.Marker = 'none';
hc(jj).Children.MarkerIndices = 3:8:numel(hc(jj).Children.MarkerIndices);
end
end
end
end
end