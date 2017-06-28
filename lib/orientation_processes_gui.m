function [L_proc,vhm_vert] = orientation_processes_gui(L_proc, L_body);
    %find the further point on x and consider the point on L_body with
    %y_min and consider the line passing these two points. Evaluate than
    %the angle 
    
    am_mean_body = mean(L_body);
    body_vh = vhm_vert(vhm_vert(:,2)>0,:);
    vh_mean_body = mean(body_vh);
    %move the vh vert
    dist = vh_mean_body - am_mean_body;
    dist = [0 dist(2) 0];
    vhm_vert = vhm_vert - dist;
    %plot_matrix(vhm_vert,'.y');
    
    %plot_matrix(L_body,'.b');
    
    %----------     open datasets
    %visual human to do the comparison
    %fileName_vhm = 'C:/Users/lavecchc/Dropbox/Internship_Melbourne/accuracy_vh_Carol/L1_vh2.stl';
    %vhm = import_STL_txt(fileName_vhm); %visual human model
    %vhm_vert = vhm.solidVertices{1,1};
    %plot_matrix(vhm_vert,'.r'); plot_matrix(L_body,'.r');   %function to plot matrix
    
    % processes divided in left and right to help in the rotations
    L_proc_dx = L_proc(L_proc(:,1)>0,:);
    L_proc_sx = L_proc(L_proc(:,1)<0,:);
    
    %------------       select points on the vh: vh_body
    vh_body = vhm_vert(vhm_vert(:,2)>-20,:);
    eps = 1;
    vh_xmean = vh_body(vh_body(:,1)<mean(vh_body(:,1))+eps & vh_body(:,1)>mean(vh_body(:,1))-eps,:);
    %plot3(vh_xmean(:,1),vh_xmean(:,2),vh_xmean(:,3),'Ob'),hold on;
    %min point on the vertebral body
    [vh_min_body_y,ind] = min(vh_xmean(:,2));
    vh_min_body = vh_xmean(ind,:);      
    
    %select the furthest point on the lateral process
    [vh_min_x,ind] = min(vhm_vert(:,1));
    vh_min_proc = vhm_vert(ind,:);
    %plot_matrix(vh_min_proc,'Ob');
     %plot_matrix(vhm_vert,'.r');

    %-------------       select points on the average model (am)
    eps = 1;
    am_xmean = L_body(L_body(:,1)<mean(L_body(:,1))+eps & L_body(:,1)>mean(L_body(:,1))-eps,:);
    %plot_matrix(am_xmean,'Ob');
    %min point on the vertebral body
    [am_min_body_y,ind] = min(am_xmean(:,2));
    am_min_body = am_xmean(ind,:);      
    
    %select the furthest point on the lateral process
    [am_min_x,ind] = min(L_proc_sx(:,1));
    am_min_proc = L_proc_sx(ind,:);
    
    %plot_matrix(am_min_proc,'Ob');
     %plot_matrix(L_proc_dx,'.r');

    %evaluate the distance between the processes 
    dist_process = am_min_proc - vh_min_proc;
    L_proc_dx = L_proc_dx - dist_process;
    L_proc_sx = L_proc_sx - dist_process;
    
    %update am_min_proc
    am_min_proc = am_min_proc - dist_process;
    %   plot_matrix(L_proc_dx,'.b');
    
    
    %------------- find vectors to evaluate the angles
    vh_bod_proc = [vh_min_body-vh_min_proc]; %vector between body and lat processes
    vh_y = [vh_min_body(:,1)-vh_min_proc(:,1) 0 vh_bod_proc(3)];
    
    am_bod_proc = [am_min_body-am_min_proc]; %vector between body and lat processes
    am_y = [am_min_body(:,1)-am_min_proc(:,1) 0 am_bod_proc(3)];
    %  dot product
    %angl_vh = acos(vh_y.*vh_bod_proc/dot(vh_y*vh_bod_proc))
    
    %   angl_am = acos(am_y.*am_bod_proc/dot(am_y*am_bod_proc))
    ang_am = atan2(norm(cross(am_y,am_bod_proc)),dot(am_y,am_bod_proc));
    ang_vh = atan2(norm(cross(vh_y,vh_bod_proc)),dot(vh_y,vh_bod_proc));
    
    %evaluate the difference between the angles and apply to a rotation
    %matrix
    ang = ang_vh - ang_am;
    ang_diff = ang*180/pi;
    
    %plot_matrix(L_proc_dx,'.b')
    %plot_matrix(L_proc_sx,'.b')

    
    %before to rotatr, translate the processes in the origin
    [L_proc_dx, L_proc_sx] = move_rotate_processes(ang, L_proc_dx, L_proc_sx);
    %plot_matrix(L_proc_dx,'.g')
    %plot_matrix(L_proc_sx,'.g')
    L_proc = [L_proc_dx; L_proc_sx];
    
    %plot_matrix(L_proc,'ob')

    
    
end
