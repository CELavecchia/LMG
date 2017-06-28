%function called from stl_test_nov.m
function [L_proc_dx,L_proc_sx,L_ped_dx, L_ped_sx] = divide_dx_sx(L_proc, L_ped)
%function [L_proc_dx,L_proc_sx, L_ped_dx, L_ped_sx] = divide_dx_sx(L_proc, L_ped)

%save in 2 different files the left and right sides
%figure, plot3(L_proc_sx(:,1),L_proc_sx(:,2),L_proc_sx(:,3),'.r'), grid on, xlabel('x'),ylabel('y'),zlabel('z'),hold on;

% to close the inferior part
L_proc_dx = L_proc(find(L_proc(:,1)<0),:);
L_proc_sx = L_proc(find(L_proc(:,1)>0),:);
%{
point_grid_dx = processing_proc(L_proc_dx);
L_proc_dx = [L_proc_dx; point_grid_dx];
%}
%point_grid_sx = processing_proc(L_proc_sx);
%point_grid_sup_sx = processing_proc_sup(L_proc_sx);
%L_proc_sx = [L_proc_sx; point_grid_sx; point_grid_sup_sx];

count = 1;
for( j = 1: length(L_ped))
   if( L_ped(j,1)>0 )
       L_ped_sx(count,:) = [L_ped(j,1),L_ped(j,2),L_ped(j,3)];
       count = count +1;
   end
end
count = 1;
for( j = 1: length(L_ped))
   if( L_ped(j,1)<0 )
       L_ped_dx(count,:) = [L_ped(j,1),L_ped(j,2),L_ped(j,3)];
       count = count +1;
   end
end

end