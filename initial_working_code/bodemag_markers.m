%% Original
rng(200,'twister')
sys1 = rss(8,2,2);
[sys2,~] = balred(sys1,2,'StateProjection','truncate');
sys3 = sys2*makeweight(1,[1e3,db2mag(-10)],db2mag(-40));
bodemag(sys1,'k',sys2,'-.r',sys3,'--b',{1e1,1e5})

%% Original Answer
rng(200,'twister')
sys1 = rss(8,2,2);
[sys2,~] = balred(sys1,2,'StateProjection','truncate');
sys3 = sys2*makeweight(1,[1e3,db2mag(-10)],db2mag(-40));
bodemag(sys1,'k',sys2,'-.r',sys3,'--b',{1e1,1e5})
hax = get(gcf,'Children');
for ii = 1:numel(hax)
    if strcmp(hax(ii).Type,'axes')
        hc = hax(ii).Children;
        for jj = 1:numel(hc)
            if strcmp(hc(jj).Type,'hggroup')
                % all Line poperties accessible here
                % set Marker and space ever 10th data point
                hc(jj).Children.Marker = 'o';
                hc(jj).Children.MarkerIndices = 1:10:numel(hc(jj).Children.MarkerIndices);
            end
        end
    end
end

%% Modified answer
rng(200,'twister')
sys1 = rss(8,2,2);
[sys2,~] = balred(sys1,2,'StateProjection','truncate');
sys3 = sys2*makeweight(1,[1e3,db2mag(-10)],db2mag(-40));
bodemag(sys1,'k',sys2,'-.r',sys3,'--b',{1e1,1e5})
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
                    hc(jj).Children.MarkerIndices = 2:11:numel(hc(jj).Children.MarkerIndices);
                end
                if jj == 3
                    % all Line poperties accessible here
                    % set Marker and space ever 10th data point
                    hc(jj).Children.Marker = 'd';
                    hc(jj).Children.MarkerIndices = 3:12:numel(hc(jj).Children.MarkerIndices);
                end
                if jj == 4
                    % all Line poperties accessible here
                    % set Marker and space ever 10th data point
                    hc(jj).Children.Marker = 'none';
                    hc(jj).Children.MarkerIndices = 3:12:numel(hc(jj).Children.MarkerIndices);
                end
            end
        end
    end
end