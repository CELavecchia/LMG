% move the bodies on the curve
%function [CM, CM_IVD,mesh_struct_IVD2,L_no_r] = position_curve_gui2(L_no_r, mesh_struct_IVD2, CM,CM_IVD, alpha, hL,rot_angle, IVD, heighV)

function [CM,CM_IVD, mesh_struct_IVD2, L_no_r, EP] = position_curve_gui2(L_no_r, mesh_struct_IVD2,  CM,CM_IVD, alpha, hL,rot_angle, IVD, heighV, EP)
  

alpha_lum = alpha;;
%to obtain more points

alpha_lum_r = pi*alpha_lum/180; %angle in radians

% total length
c = sum(hL) + sum(IVD); %(mm)
%radius corresponding circumference
r = c /( 2*sin(alpha_lum) );

eps = 0.6;
y_cent = 0;
z_cent = 0;
th = 0: alpha_lum_r/100: alpha_lum_r;


th = 0: (pi/3)/1000 : pi/3;
yunit = -r * cos(th) + y_cent;
zunit = r * sin(th) + z_cent;
x = zeros(length(zunit)); 

%{
for(j=1:1)
  plot3(mesh_struct_IVD2(j).V(:,1),mesh_struct_IVD2(j).V(:,2), mesh_struct_IVD2(j).V(:,3),'.b'),hold on;
  plot3(CM_IVD(:,1), CM_IVD(:,2), CM_IVD(:,3),'Or'),hold on;
  
end
%}
%plot3(CM(1,1),CM(1,2),CM(1,3),'Og'),hold on;

%figure(2)

%plot3(x,yunit, zunit, '*b'),grid on, xlabel('x'),hold on, ylabel('y'), zlabel('z');


% evaluate the height and its half, evaluate the same distance on the curve
% and find the point
half_h = heighV./2;
%rot_angle
%r1 = half_h(1) /( 2*sin(rot_angle(1)) );

%rot_angle is the angle of the cross product!!! 
%y_gc = -r*cos(rot_angle(1)) + y_cent
%z_gc = r*sin(rot_angle(1)) + y_cent
count_v =1;
z1=zunit(1);
z_ip = 0; count_ang = 1;
for(body=1:5)
   % plot3(L_no_r(count_v).VT(:,1),L_no_r(count_v).VT(:,2),L_no_r(count_v).VT(:,3),'.b'),hold on;
   % plot3(CM(body,1),CM(body,2),CM(body,3),'Og'),hold on;

        %if(mod(body,2) ~= 0) %if vert body
  
         if(count_v==1)
     
            y_ip = -r ;
            %move vertebrae; traslation vector
            t(count_v,:) = [ y_ip-CM(count_v,8) 0] ; %keep the same z of CM

            L_no_r(count_v).VT = [L_no_r(count_v).VT(:,1) L_no_r(count_v).VT(:,2)+(t(count_v,1)) L_no_r(count_v).VT(:,3)+5]; 
            EP(count_v).sup = [EP(count_v).sup(:,1) EP(count_v).sup(:,2)+(t(count_v,1)) EP(count_v).sup(:,3)+5];
            EP(count_v).inf = [EP(count_v).inf(:,1) EP(count_v).inf(:,2)+(t(count_v,1)) EP(count_v).inf(:,3)+5];
            
            %update the CM
            CM(count_v,2:3) = [ CM(count_v,2)+t(count_v,1) CM(count_v,3)]; % z in the first one stay the same
            CM(count_v,5:6) = [ CM(count_v,5)+t(count_v,1) CM(count_v,6)]; % z in the first one stay the same
            CM(count_v,8:9) = [ CM(count_v,8)+t(count_v,1) CM(count_v,9)]; % z in the first one stay the same
            z_old = CM(count_v,9);

            z1 = CM(1,9) - IVD(count_v) - half_h(count_v);
            %plot3(L_no_r(count_v).VT(:,1),L_no_r(count_v).VT(:,2),L_no_r(count_v).VT(:,3),'.r'),hold on;
            count_v = count_v +1;
            %z1 = CM(1,9);    
        
         elseif count_v==2
            z1 = z1-half_h(count_v);
            z_ip = z_old-(z1);

            angle_cor(count_ang) = atan(z_ip/r);

            z1_new = r*sin(angle_cor(count_ang))+abs(z_old);
            z_ip = z1_new;
            %abs(r)
            %r*cos(angle_cor)
            y1_new = r*cos(angle_cor(count_ang));

            y_ip = -y1_new;% -r * cos(angle_cor) + y_cent
            %move vertebrae; traslation vector
                %CM
                
            t(count_v,:) = [ y_ip-abs(CM(count_v,8)) z_ip-abs(CM(count_v,9))];

            L_no_r(count_v).VT = [L_no_r(count_v).VT(:,1) L_no_r(count_v).VT(:,2)+(t(count_v,1)) L_no_r(count_v).VT(:,3)-(t(count_v,2))];
            EP(count_v).sup = [EP(count_v).sup(:,1) EP(count_v).sup(:,2)+(t(count_v,1)) EP(count_v).sup(:,3)-(t(count_v,2))];
            EP(count_v).inf = [EP(count_v).inf(:,1) EP(count_v).inf(:,2)+(t(count_v,1)) EP(count_v).inf(:,3)-(t(count_v,2))];
            
            %update CM
            CM(count_v,2:3) = [ CM(count_v,2)+t(count_v,1) CM(count_v,3)-t(count_v,2)]; % z in the first one stay the same
            CM(count_v,5:6) = [ CM(count_v,5)+t(count_v,1) CM(count_v,6)-t(count_v,2)]; % z in the first one stay the same
            CM(count_v,8:9) = [ CM(count_v,8)+t(count_v,1) CM(count_v,9)-t(count_v,2)]; % z in the first one stay the same
                
            
           z_old = z_ip;

            z1 = z1 - IVD(count_v) - half_h(count_v); %z1 evaluated in this loop + the IVD height + the half body of this loop

            count_ang = count_ang +1;
            %plot3(L_no_r(count_v).VT(:,1),L_no_r(count_v).VT(:,2),L_no_r(count_v).VT(:,3),'.r'),hold on;
            count_v = count_v +1;
    else
            z1 = z1 - half_h(count_v);
            z = abs(z1) - z_old; %length segment
            z_perp = z*cos(angle_cor(count_ang-1));     %perpendicular segment

            angle_cor(count_ang) = atan(z_perp/r);

            z_ip = z_ip + r*sin(angle_cor(count_ang));
            y_ip = -r*cos(sum(angle_cor));%(count_ang)+angle_cor(count_ang-1) )

            t(count_v,:) = [ y_ip-abs(CM(count_v,8)) z_ip-abs(CM(count_v,9))];
            L_no_r(count_v).VT = [L_no_r(count_v).VT(:,1) L_no_r(count_v).VT(:,2)+(t(count_v,1)) L_no_r(count_v).VT(:,3)-(t(count_v,2))];
            EP(count_v).sup = [EP(count_v).sup(:,1) EP(count_v).sup(:,2)+(t(count_v,1)) EP(count_v).sup(:,3)-(t(count_v,2))];
            EP(count_v).inf = [EP(count_v).inf(:,1) EP(count_v).inf(:,2)+(t(count_v,1)) EP(count_v).inf(:,3)-(t(count_v,2))];
            
            
            
             CM(count_v,2:3) = [ CM(count_v,2)+t(count_v,1) CM(count_v,3)-t(count_v,2)]; % z in the first one stay the same
            CM(count_v,5:6) = [ CM(count_v,5)+t(count_v,1) CM(count_v,6)-t(count_v,2)]; % z in the first one stay the same

            CM(count_v,8:9) = [ CM(count_v,8)+t(count_v,1) CM(count_v,9)-t(count_v,2)]; % z in the first one stay the same
          
            
           % plot3(L_no_r(count_v).VT(:,1),L_no_r(count_v).VT(:,2),L_no_r(count_v).VT(:,3),'.r'),hold on;

            z_old = z_ip;
            if count_v ~= 5
                z1 = z1 - IVD(count_v) - half_h(count_v);
            end
                count_ang = count_ang +1;
            %t(j,:) = [ y_ip-abs(CM(j,8)) z_ip-abs(CM(j,9))];
            count_v = count_v +1;
         
         end
  
         
end




%move IVD
for(j =1:4)

  t = [ CM(j,1)-CM_IVD(j,4) CM(j,2)-CM_IVD(j,5) CM(j,3)-CM_IVD(j,6)];
  
  %plot3(mesh_struct_IVD2(j).V(:,1),mesh_struct_IVD2(j).V(:,2), mesh_struct_IVD2(j).V(:,3),'.g'),hold on;
   CM_IVD(j,1:3) =[CM_IVD(j,1)-t(1),CM_IVD(j,2)+abs(t(2)),CM_IVD(j,3)-abs(t(3))];
   CM_IVD(j,4:6) =[CM_IVD(j,4)-t(1),CM_IVD(j,5)+abs(t(2)),CM_IVD(j,6)-abs(t(3))];
   CM_IVD(j,7:9) =[CM_IVD(j,7)-t(1),CM_IVD(j,8)+abs(t(2)),CM_IVD(j,9)-abs(t(3))];

   mesh_struct_IVD2(j).V = [ mesh_struct_IVD2(j).V(:,1)-t(1) mesh_struct_IVD2(j).V(:,2)+abs(t(2)) mesh_struct_IVD2(j).V(:,3)-abs(t(3))];
   
      
   
  %plot3(CM_IVD(j,1),CM_IVD(j,2),CM_IVD(j,3),'.c'),hold on;
  %plot3(mesh_struct_IVD2(j).V(:,1),mesh_struct_IVD2(j).V(:,2), mesh_struct_IVD2(j).V(:,3),'.g'),hold on;
end



[CM, CM_IVD,mesh_struct_IVD2,L_no_r, EP] = check_position_arch_gui2(CM, CM_IVD, mesh_struct_IVD2, L_no_r, EP);
%[CM,mesh_struct_IVD2,L_no_r, EP] = check_position_arch_gui2(CM, mesh_struct_IVD2, L_no_r, EP);



end







