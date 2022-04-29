%% bc1.m
% For systems 3,11,19 (second, sixth and tenth (last) beams)
for j = 1:ns
    if  j == 3 || j == 11 || j == 19
        if mod(j,2) == 1
            figure()
            bodemag(G(j).W1,'k',G(j).W1e,'-.r',G(j).W1s,'--g', {1e1,1e5})
        end
    end
end
% export_fig weights_c1_3.svg -painters
% export_fig weights_c1_11.svg -painters
% export_fig weights_c1_19.svg -painters


%% bc2.m
% For systems 3,7,13 (second, fourth and seventh (last) beams)
for j = 1:ns
    if  j == 3 || j == 7 || j == 13
        if mod(j,2) == 1
            figure()
            bodemag(G(j).W1,'k',G(j).W1e,'-.r',G(j).W1s,'--g', {1e1,1e5})
        end
    end
end
% export_fig weights_c2_3.svg -painters
% export_fig weights_c2_7.svg -painters
% export_fig weights_c2_13.svg -painters

%% bc3.m
% For systems 3,5,9 (second, third and fifth (last) beams)
for j = 1:ns
    if  j == 3 || j == 5 || j == 9
        if mod(j,2) == 1
            figure()
            bodemag(G(j).W1,'k',G(j).W1e,'-.r',G(j).W1s,'--g', {1e1,1e5})
        end
    end
end
% export_fig weights_c3_3.svg -painters
% export_fig weights_c3_5.svg -painters
% export_fig weights_c3_9.svg -painters
