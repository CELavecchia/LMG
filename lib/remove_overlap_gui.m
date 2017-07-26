%{
    remove_overlap: function to remove doubled points
    %}

function [L_body, L_proc, L_ped, lam_new, EPsup,EPinf] =...
   remove_overlap_gui(L_body, L_proc_grid_nt, L_lam_nt,EPsup,EPinf,...
    EPWu_half, EPDu, hL, lam_l, sc_d, sc_w, PDW, PDH, TP_wu, TP_wi,...
   PDt,PDs);


%first at all mirror it
L_body= [L_body(:,1), L_body(:,2), L_body(:,3); ...
    -L_body(:,1), L_body(:,2), L_body(:,3)];
L_lam_nt = [L_lam_nt(:,1), L_lam_nt(:,2), L_lam_nt(:,3); ...
    -L_lam_nt(:,1), L_lam_nt(:,2), L_lam_nt(:,3)];



[L_body,L_proc, L_ped, L_lam, EPsup,EPinf] = position_posterior_el_gui(L_body,L_lam_nt,L_proc_grid_nt,EPsup,EPinf,...
    EPWu_half, EPDu, hL,lam_l, sc_d, sc_w, PDW, PDH, TP_wu, TP_wi, PDt, PDs );
%{
[L_body,L_proc, L_ped, L_lam] = position_posterior_el(L_body,L_lam_nt,L_proc_grid_nt,L_ped_nt,...
    EPWu_half, EPDu, hL,lam_l, sc_d, sc_w, PDW, PDH, TP_wu, TP_wi);%, vh_vert
%}
%plot3(L_ped(:,1), L_ped(:,2), L_ped(:,3),'.b'),hold on;
%plot3(L_body(:,1), L_body(:,2), L_body(:,3),'.r'),hold on;

%divide in two different dataasedts, left (sx) and rights(dx) ones
[L_proc_dx,L_proc_sx, L_ped_dx, L_ped_sx] = divide_dx_sx_gui(L_proc, L_ped);



ped_dx_new = check_overlap_gui( L_body,  L_ped_dx );
ped_sx_new = check_overlap_gui( L_body,  L_ped_sx );

ped_dx_new = check_overlap_ped_proc_gui( L_proc_dx, ped_dx_new  );   
ped_sx_new = check_overlap_ped_proc_gui( L_proc_sx, ped_sx_new  );


lam_new = check_overlap_lam_gui( L_proc_dx,  L_lam );   
lam_new = check_overlap_lam_gui( L_proc_sx,  lam_new );   

%put again the proc and ped datasets together:
L_ped = [ped_dx_new; ped_sx_new];


%plot3(L_ped(:,1), L_ped(:,2), L_ped(:,3),'.b'),hold on;
%plot3(L_body(:,1), L_body(:,2), L_body(:,3),'.r'),hold on;


end
