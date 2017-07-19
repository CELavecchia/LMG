% function to positioning vertebrae and IVD in the right positions using
% the CM point
%clc; clear; close all;
%function  [vert_t, L_lam_t, L_proc_t, L_ped_t, IVD_data_t] = position_curve(vert, L_lam, L_proc, L_ped, IVD_data, CM, CM_IVD)
function  [L_no, mesh_struct_IVD2] = position_curve_gui(L_no_r, mesh_struct_IVD2, CM, alpha, hL, IVD)
  
%CM_IVD = 

%circumference
alpha_lum = alpha;%43.49;
%to obtain more points
%alpha_lum = 45;
alpha_lum_r = pi*alpha_lum/180; %angle in radians

% total length
c = sum(hL) + sum(IVD); %(mm)
%radius corresponding circumference
r = c /( 2*sin(alpha_lum) );

%
y_cent = 0;
z_cent = 0;
th = 0: alpha_lum_r/100: alpha_lum_r;


th = 0: (pi/3)/100 : pi/3;
yunit = -r * cos(th) + y_cent;
zunit = r * sin(th) + z_cent;
x = zeros(length(zunit)); 
plot3(x,yunit, zunit, '*b'),grid on, xlabel('x'), ylabel('y'), zlabel('z');
%hold on;
%CM
% first point -> yunit(1),zunit(1);
%height vertebrae and IVD, diff sup and inf points only along z
dist_v = [] ; dist_i = [];
for(i = 1:5)
   
    dist_v(i) = sqrt( (CM(i,6)-CM(i,3))^2 );
    %dist_i(i) = sqrt( (CM_IVD(i,6)-CM_IVD(i,3))^2 );
    
end

% find the point on the curve distant from the first one of the length
% dist_v
eps = [0.1*dist_v; 0.1*dist_i];
eps1 = 1;
eps2 = 2;
y_cm_c = []; z_cm_c = [];
y_pos(1) = yunit(1);
z_pos(1) = zunit(1);
count = 1; count_v = 1; count_i = 1;
for(body = 1:9)        
    d=0;
    if(mod(body,2) ~= 0)
        %vert
        for(i = 1:length(zunit))
            d = sqrt( (y_pos(count)-yunit(i))^2 + (z_pos(count)-zunit(i))^2 );
            if ( d>= dist_v(count_v)-eps1 && d<= dist_v(count_v)+eps1 ) %find the point corrispondent to that area
                y_cm_c = yunit(i);
                z_cm_c = zunit(i);
            end
        end
        %plot3(zeros(size(y_cm_c)),y_cm_c,z_cm_c, 'Ob'), hold on;
        count_v = count_v +1;
    end
    %{
    if(mod(body,2) == 0)
        %ivd
        for(i = 1:length(zunit))
            d = sqrt( (y_pos(count)-yunit(i))^2 + (z_pos(count)-zunit(i))^2 );
            if ( d>= dist_i(count_i)-eps1 && d<= dist_i(count_i)+eps1 ) %find the point corrispondent to that area
                y_cm_c = yunit(i);
                z_cm_c = zunit(i);
            end
        end
       % plot3(zeros(size(y_cm_c)),y_cm_c,z_cm_c, 'Or')
       count_i = count_i +1;
    end
    %}
    % in this way I found the points in corrispondence of the end of each
    % body in the case of perfect contact
    count = count +1; % it s right to put first the counter, since the first place is already taken!
    y_pos(count) = y_cm_c;%+eps2;
    z_pos(count) = z_cm_c;%-eps2;
    % giving space
    
end
%plot(y_pos, z_pos,'^b')
size(y_pos); % 10 elements-> ok.

%size(yunit)

% y_pos and z_pos are the positions of the ending and beginning of each
% body. now I need the middle point and I need to move them a little
count = 2;
d = 0; 
eps = 0;
eps_seek = 2;
CM_y = []; CM_z = [];
index_inf = (1:(length(y_pos)-1));
index_sup = (2:(length(y_pos)));

for(i = 1:(length(y_pos)-1))
      %distance between inf and sup points 
   d_points(i) = sqrt( (y_pos(index_sup(i))-y_pos(index_inf(i)))^2 + (z_pos(index_sup(i))-z_pos(index_inf(i)))^2 ) ;
 
end
count_dis = 1;


for(i = 1:length(y_pos)-1)
   %evaluate the distance between the points previously saved
   for(j = 1 : length(yunit) )  % evaluate the points on the curve, and find the one at the distance
       d = sqrt( (y_pos(i)-yunit(j))^2 + (z_pos(i)-zunit(j))^2 );
         %d2 = sqrt( (CM_y(i)-yunit(j))^2 + (CM_z(i)-zunit(j))^2 );
       %d = sqrt( (y_tmp-yunit(j))^2 + ( z_tmp - zunit(j) )^2 );
            if ( d>= (d_points(count_dis)+eps)/2-eps_seek && ...
                    d<= (d_points(count_dis)+eps)/2+eps_seek && zunit(j)<z_pos(i))% && ...
                  % d2>= (d_points(count_dis)+eps)/2-eps_seek && d2<= (d_points(count_dis)+eps)/2+eps_seek ) %find the point corrispondent to that area
            %if ( d>= (d_points(i)+eps)/2-eps_seek && d<= (d_points(i)+eps)/2+eps_seek && zunit(j)<z_pos(i)) %find the point corrispondent to that area
    
                CM_y(i) = yunit(j);
                CM_z(i) = zunit(j);
                  %CM_y(i+1) = yunit(j);
                %CM_z(i+1) = zunit(j);
            end
   end
   %plot( CM_y, CM_z, '+m'), grid on;
   %y_tmp = CM_y(i+1) ;
   %z_tmp = CM_z(i+1) ;
   eps = eps + 8;
   count = count + 1;
   count_dis = count_dis+1;
end
l = length(CM_y);

x = zeros(length(CM_y));


%hold on;
%size(vert)
%calculate the traslation vector between the pointsCM_y CM_z and the CM of
%vertebrae and IVD
count_v = 1;
count_i = 1;
count_bodv = 1;
count_bodi = 1;
vert_t = [] ;
IVD_data_t = [];
count_surf = 1;
L_lam_t = []; L_proc_t = []; L_ped_t = [];
for(body = 1:9)
   
    if(mod(body,2) ~= 0) %it s a vertebrae
        t = [CM(count_bodv,8)-CM_y(body), CM(count_bodv,9)-CM_z(body)];
        
        L_no(count_bodv).VT = [ L_no_r(count_bodv).VT(:,1), L_no_r(count_bodv).VT(:,2)-t(1), L_no_r(count_bodv).VT(:,3)-t(2)];
                
        
                  
        count_v = count_v +3;
        count_bodv = count_bodv + 1 ;
        count_surf = count_surf+1;
        %{
    elseif (mod(body,2) == 0)
        t = [CM_IVD(count_bodi,8)-CM_y(body), CM_IVD(count_bodi,9)-CM_z(body)];
  
        %mesh 
         anulus(count_i).vertices =  [anulus(count_i).vertices(:,1), ...
            anulus(count_i).vertices(:,2)-t(1), anulus(count_i).vertices(:,3)-t(2)] ;
         nucleus(count_i).vertices =  [ nucleus(count_i).vertices(:,1), ...
            nucleus(count_i).vertices(:,2)-t(1), nucleus(count_i).vertices(:,3)-t(2)] ;
      
        count_i = count_i +1;
        count_bodi = count_bodi +1;
        %}
    end

    
end





