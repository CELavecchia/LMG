function L_ped =orient_pedicles_gui(L_ped)
   % plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.b'),hold on;

   p = polyval(L_ped(:,2),L_ped(:,3),1);
   mean_pos = mean(L_ped);
   %traslate in zero
   L_ped_t = [L_ped(:,1)-mean_pos(:,1), L_ped(:,2)-mean_pos(:,2), L_ped(:,3)-mean_pos(:,3)];
   %plot3(L_ped_t(:,1),L_ped_t(:,2),L_ped_t(:,3),'.r'),hold on;
   
   angl = -20*pi/180; %value that can be imported
   %rotate around x
   rot_x = [1 0 0; 0 cos(angl) sin(angl); 0 -sin(angl) cos(angl)];
   L_ped_r = L_ped_t*rot_x;
   %  plot3(L_ped_r(:,1),L_ped_r(:,2),L_ped_r(:,3),'.b'),hold on;
    %move it back
   
    L_ped_dx = L_ped_r(L_ped_r(:,1)<0,:); mean_p_dx = mean(L_ped_dx);
    L_ped_dx = [L_ped_dx(:,1)-mean_p_dx(:,1),L_ped_dx(:,2:3)];
       %plot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'.g')

    L_ped_sx = L_ped_r(L_ped_r(:,1)>0,:); mean_p_sx = mean(L_ped_sx);
    L_ped_sx = [L_ped_sx(:,1)-mean_p_sx(:,1),L_ped_sx(:,2:3)];
   % plot3(L_ped_sx(:,1),L_ped_sx(:,2),L_ped_sx(:,3),'.g')
    %in z
    angl = -10*pi/180; %this angle could be imported as well
    rot_z = [ cos(angl) sin(angl) 0; -sin(angl) cos(angl) 0; 0 0 1];
      L_ped_sx = L_ped_sx*rot_z;
      angl = 10*pi/180; %this angle could be imported as well
    rot_z = [ cos(angl) sin(angl) 0; -sin(angl) cos(angl) 0; 0 0 1];
      L_ped_dx = L_ped_dx*rot_z;

     L_ped_dx = [L_ped_dx(:,1)+mean_p_dx(:,1),L_ped_dx(:,2:3)];
     L_ped_sx = [L_ped_sx(:,1)+mean_p_sx(:,1),L_ped_sx(:,2:3)];
    L_ped = [L_ped_dx;L_ped_sx];
   %plot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'.m')
   %plot3(L_ped_r(:,1),L_ped_r(:,2),L_ped_r(:,3),'.r'),hold on;
    %move it back
    L_ped = [L_ped(:,1)+mean_pos(:,1), L_ped(:,2)+mean_pos(:,2), L_ped(:,3)+mean_pos(:,3)];
    
  %  plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.r'),hold on;

    
end