function [L_body, L_ped, L_proc, L_lam2] = double_check_positions_el_gui(L_body, L_ped, L_proc, L_lam2, sc_d, sc_w, PD_w, PD_h);

% remember: if you move the lamina, move all together with the processes

%1 sc_d--- lamina 
evaluate_body = L_body(L_body(:,1)>=-0.05 & L_body(:,1)<=0.05,: );
[min_body, ind] = min(evaluate_body(:,2));

evaluate_lam = L_lam2(L_lam2(:,1)>=-0.05 & L_lam2(:,1)<=0.05,: );
[max_lam, ind] = max(evaluate_lam(:,2));

dist_init = abs(max_lam) - abs(min_body);
t_y = sc_d - dist_init; % traslate in order to obtain sc_d

%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.b'),hold on;
%plot3(L_lam2(:,1),L_lam2(:,2),L_lam2(:,3),'.b'),hold on;

L_lam2 = [L_lam2(:,1), L_lam2(:,2)-t_y, L_lam2(:,3)];
L_proc = [L_proc(:,1), L_proc(:,2), L_proc(:,3)];%-t_y
%plot3(L_lam2(:,1),L_lam2(:,2),L_lam2(:,3),'.r'),hold on;

% sc_w ---> pedicles
eps = 1;
L_ped_mean = L_ped(L_ped(:,1)<0 & L_ped(:,2)<=mean(L_ped(:,2))+eps & L_ped(:,2)>=mean(L_ped(:,2))-eps,:); % I consider only one of the pedicles
mean_x = mean(L_ped_mean(:,1));
%plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.b'),hold on;

scw_init = abs(mean_x) - PD_w/2;%+2;  
t_x = (sc_w/2) - scw_init;
L_ped_min = L_ped(L_ped(:,1)<0,:);
L_ped_min = [L_ped_min(:,1)-t_x, L_ped_min(:,2), L_ped_min(:,3)];
L_ped_max = L_ped(L_ped(:,1)>0,:);
L_ped_max = [L_ped_max(:,1)+t_x, L_ped_max(:,2), L_ped_max(:,3)];

L_ped = [L_ped_min; L_ped_max];
%plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.r'),hold on;


%sc_d ----> pedicles
mean_y = mean(L_ped_mean(:,2)); %I want the mean of the pedicles to be in corrispondence of the mean of the sc_d
scd_dist = abs(mean_y)- abs(min_body-sc_d/2);
L_ped = [L_ped(:,1), L_ped(:,2)-scd_dist, L_ped(:,3)];
%plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.b'),hold on;

end