%rotate the IVD 
function [vert_r , CM ] = moving_rot_post_mesh_ivd_gui(vert, CM, t, mat)    
  
         vert_m(:,1:3) = [vert(:,1), vert(:,2)+t(1), vert(:,3)+abs(t(2))];
  
        
        CMt(8:9) = [CM(8)+t(1), CM(9)+abs(t(2))];
        CMt(5:6) = [CM(5)+t(1), CM(6)+abs(t(2))];
        CMt(2:3) = [CM(2)+t(1), CM(3)+abs(t(2))];
        
        vert_r(:,2:3) = vert_m(:,2:3)*mat;
      
                  
        % ok plot3(vert_m(:,count_colv),vert_r(:,count_colv+1), vert_r(:,count_colv+2),'.r');
        CM_r(8:9) =  CMt(8:9)*mat;
        CM_r(5:6) =  CMt(5:6)*mat;
        CM_r(2:3) = CMt(2:3)*mat;
        
        %update traslation
        CM(1:3) = [ CM(1), -CM_r(2)-t(1),-CM_r(3)-abs(t(2))];
        CM(4:6) = [ CM(4), -CM_r(5)-t(1),-CM_r(6)-abs(t(2))];
        CM(7:9) = [ CM(7), -CM_r(8)-t(1),-CM_r(9)-abs(t(2))];

        
        vert_r(:,1:3) = [vert_m(:,1), -vert_r(:,2)-t(1), -vert_r(:,3)-abs(t(2))];
  
         
        
end