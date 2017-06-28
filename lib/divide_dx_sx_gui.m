%function called from stl_test_nov.m
function [L_proc_dx,L_proc_sx,L_ped_dx, L_ped_sx] = divide_dx_sx_gui(L_proc, L_ped)
%function [L_proc_dx,L_proc_sx, L_ped_dx, L_ped_sx] = divide_dx_sx(L_proc, L_ped)

%save in 2 different files the left and right sides
%figure, plot3(L_proc_sx(:,1),L_proc_sx(:,2),L_proc_sx(:,3),'.r'), grid on, xlabel('x'),ylabel('y'),zlabel('z'),hold on;

% to close the inferior part
L_proc_dx = L_proc(L_proc(:,1)<0,:);
L_proc_sx = L_proc(L_proc(:,1)>0,:);


L_ped_sx = L_ped(L_ped(:,1)>0,:);
L_ped_dx = L_ped(L_ped(:,1)<0,:);


end