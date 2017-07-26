%function [vert_r , CM ] = moving_rot(vert, CM, count_colv, count_v, t, mat)%, CM_y, CM_z)    
 function [vert_r ,EPsup_r, EPinf_r,  CM ] = moving_rot_post_mesh_gui(vert, CM, EPsup, EPinf, t, mat)    
  
 


        vert_m(:,1:3) = [vert(:,1), vert(:,2)+t(1), vert(:,3)+abs(t(2))];
        EPsup(:,1:3) = [EPsup(:,1), EPsup(:,2)+t(1), EPsup(:,3)+abs(t(2))];
        EPinf(:,1:3) = [EPinf(:,1), EPinf(:,2)+t(1), EPinf(:,3)+abs(t(2))];

        %plot3(vert(:,1),vert(:,2),vert(:,3),'.r'),hold on
        %plot3(CM(1,7),CM(1,8), CM(1,9),'oy'),hold on;
        %plot3(vert_m(:,count_colv),vert_m(:,count_colv+1), vert_m(:,count_colv+2),'.g');
        CMt(8:9) = [CM(8)+t(1), CM(9)+abs(t(2))];
        CMt(5:6) = [CM(5)+t(1), CM(6)+abs(t(2))];
        CMt(2:3) = [CM(2)+t(1), CM(3)+abs(t(2))];
        
        vert_r(:,2:3) = vert_m(:,2:3)*mat;
        EPsup_r(:,2:3) = EPsup(:,2:3)*mat;
        EPinf_r(:,2:3) = EPinf(:,2:3)*mat;

                  
        % ok plot3(vert_m(:,count_colv),vert_r(:,count_colv+1), vert_r(:,count_colv+2),'.r');
        CM_r(8:9) =  CMt(8:9)*mat;
        CM_r(5:6) =  CMt(5:6)*mat;
        CM_r(2:3) = CMt(2:3)*mat;
        
        %update traslation
        CM(1:3) = [ CM(1), -CM_r(2)-t(1),-CM_r(3)-abs(t(2))];
        CM(4:6) = [ CM(4), -CM_r(5)-t(1),-CM_r(6)-abs(t(2))];
        CM(7:9) = [ CM(7), -CM_r(8)-t(1),-CM_r(9)-abs(t(2))];

        
        %plot3(CM(count_v,1),CM(count_v,2), CM(count_v,3),'or')
        %plot3(CM(count_v,4),CM(count_v,5), CM(count_v,6),'or')
        vert_r(:,1:3) = [vert_m(:,1), -vert_r(:,2)-t(1), -vert_r(:,3)-abs(t(2))];
        EPsup_r(:,1:3) = [EPsup(:,1), -EPsup_r(:,2)-t(1), -EPsup_r(:,3)-abs(t(2))];
        EPinf_r(:,1:3) = [EPinf(:,1), -EPinf_r(:,2)-t(1), -EPinf_r(:,3)-abs(t(2))];

        %plot3(CM(1,7),CM(1,8), CM(1,9),'or'),hold on;
        %plot3(vert_r(:,1),vert_r(:,2),vert_r(:,3),'.b'),hold on;
end