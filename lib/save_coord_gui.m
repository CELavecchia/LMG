%{
 function to save points on the superior and inferior surfaces. 

1st method: saving points on the coords definition
2nd method: import surfaces and perimeters

%}
function save_coord_gui(L_no);

ind = 1 ;
for(j = 1:5)
fprintf('saving the coordinates for the vertebrae L%d \n',j)

        name_file = sprintf('./output/point_cloud/L_body_%d.txt',j);
        %name_file = sprintf('./OUTPUT_proc/vh_model/points_cloud/L_body_%d.txt',j);
        write_txt_3d_gui( L_no(j).body(:,1), L_no(j).body(:,2), L_no(j).body(:,3), name_file );
        
        name_file = sprintf('./output/point_cloud/L_lam_%d.txt',j);
        %name_file = sprintf('./OUTPUT_proc/vh_model/points_cloud/L_lam_%d.txt',j);
        write_txt_3d_gui( L_no(j).lam(:,1), L_no(j).lam(:,2), L_no(j).lam(:,3), name_file );
        
        name_file = sprintf('./output/point_cloud//L_ped_%d.txt',j);
        %name_file = sprintf('./OUTPUT_proc/vh_model/points_cloud/L_ped_%d.txt',j);
        %L_ped(:,ind:ind+2) = [L_ped_t(:,ind), L_ped_t(:,ind+1), L_ped_t(:,ind+2)];
        write_txt_3d_gui( L_no(j).ped(:,1), L_no(j).ped(:,2), L_no(j).ped(:,3), name_file );

        name_file = sprintf('./output/point_cloud//L_proc_grid_%d.txt',j);
        %name_file = sprintf('./OUTPUT_proc/vh_model/points_cloud/L_proc_grid_%d.txt',j);
        %L_proc_grid(:,ind:ind+2) = [L_proc_grid_t(:,ind), L_proc_grid_t(:,ind+1), L_proc_grid_t(:,ind+2)];
        write_txt_3d_gui( L_no(j).proc(:,1), L_no(j).proc(:,2), L_no(j).proc(:,3), name_file );
        
        ind = ind+3;
end
%plot3(sup_point(:,1), sup_point(:,2), sup_point(:,3),'.r'), xlabel('x'), ylabel('y'), zlabel('z');
%plot3(inf_point(:,1), inf_point(:,2), inf_point(:,3),'Og'), xlabel('x'), ylabel('y'), zlabel('z');

end

