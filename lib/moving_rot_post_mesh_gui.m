%function [vert_r , CM ] = moving_rot(vert, CM, count_colv, count_v, t, mat)%, CM_y, CM_z)    
 function [vert_r , CM ] = moving_rot_post_mesh_gui(vert, CM, count_colv,count_v, t, mat)    
  
        vert_m(:,1:3) = [vert(:,1), vert(:,2)+t(1), vert(:,3)+abs(t(2))];
        
        %plot3(vert_m(:,count_colv),vert_m(:,count_colv+1), vert_m(:,count_colv+2),'.g');
        CMt(5:6) = [CM(5)+t(1), CM(6)+abs(t(2))];
        CMt(2:3) = [CM(2)+t(1), CM(3)+abs(t(2))];
        
        vert_r(:,2:3) = vert_m(:,2:3)*mat;
     
                  
        % ok plot3(vert_m(:,count_colv),vert_r(:,count_colv+1), vert_r(:,count_colv+2),'.r');
        CM_r(5:6) =  CMt(5:6)*mat;
        CM_r(2:3) = CMt(2:3)*mat;
        
        %update traslation
        %t = [ CM_r(count_v,8)-CM_y, CM_r(count_v,9)-CM_z ];
        CM(1:3) = [ CM(1), -CM_r(2)-t(1),-CM_r(3)-abs(t(2))];
        CM(4:6) = [ CM(4), -CM_r(5)-t(1),-CM_r(6)-abs(t(2))];
        
        %plot3(CM(count_v,1),CM(count_v,2), CM(count_v,3),'or')
        %plot3(CM(count_v,4),CM(count_v,5), CM(count_v,6),'or')
        vert_r(:,1:3) = [vert_m(:,1), -vert_r(:,2)-t(1), -vert_r(:,3)-abs(t(2))];
        
      
end