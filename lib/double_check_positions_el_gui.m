function [L_body, L_ped, L_proc, L_lam2] = double_check_positions_el_gui(L_body, L_ped, L_proc, L_lam2, sc_d);

% remember: if you move the lamina, move all together with the processes

%1 sc_d--- lamina 
evaluate_body = L_body(L_body(:,1)>=-0.05 & L_body(:,1)<=0.05,: );
[min_body, ind] = min(evaluate_body(:,2));

evaluate_lam = L_lam2(L_lam2(:,1)>=-0.05 & L_lam2(:,1)<=0.05,: );
[max_lam, ind] = max(evaluate_lam(:,2));

dist_init = abs(max_lam) - abs(min_body);
%sc_d
t_y = sc_d - dist_init; % traslate in order to obtain sc_d
%figure
%plot3(L_lam2(:,1),L_lam2(:,2),L_lam2(:,3),'.b'),grid on;


L_lam2 = [L_lam2(:,1), L_lam2(:,2)-t_y, L_lam2(:,3)];
L_proc = [L_proc(:,1), L_proc(:,2)-t_y, L_proc(:,3)];%-t_y



end