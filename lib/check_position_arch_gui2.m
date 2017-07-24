%check the alignment between the vertebrae and IVD

function [CM, CM_IVD, mesh_struct_IVD2, L_no_r ] =check_position_arch_gui2(CM, CM_IVD, mesh_struct_IVD2, L_no_r)


%align the inferior point of the IVD with the superior of the following
%vertebrae and the inferior point of the vertebrae with the superior ones
%of the IVD
%figure(3)



%figure
%plot3(L_no_r(1).VT(:,1),L_no_r(1).VT(:,2),L_no_r(1).VT(:,3),'.r'),hold on;

% -------  vert 1 IDV1        
t = [ (CM(1,1))-(CM_IVD(1,4)) (CM(1,2))-(CM_IVD(1,5)) (CM(1,3))-(CM_IVD(1,6))]; %inferior point on the vert and superior on the ivd
         
mesh_struct_IVD2(1).V =  [mesh_struct_IVD2(1).V(:,1)+abs(t(1)) ...
  mesh_struct_IVD2(1).V(:,2)-abs(t(2)) mesh_struct_IVD2(1).V(:,3)+abs(t(3)) ];
 
CM_IVD(1,1:3) = [CM_IVD(1,1)+abs(t(1)) CM_IVD(1,2)-abs(t(2)) CM_IVD(1,3)+abs(t(3))];
CM_IVD(1,4:6) = [CM_IVD(1,4)+abs(t(1)) CM_IVD(1,5)-abs(t(2)) CM_IVD(1,6)+abs(t(3))];
CM_IVD(1,7:9) = [CM_IVD(1,7)+abs(t(1)) CM_IVD(1,8)-abs(t(2)) CM_IVD(1,9)+abs(t(3))];


% plot3(mesh_struct_IVD2(1).V(:,1),mesh_struct_IVD2(1).V(:,2),mesh_struct_IVD2(1).V(:,3),'.b'),hold on;

% -------- IVD1 - vert2

 t = [ (CM_IVD(1,1))-(CM(2,4)) (CM_IVD(1,2))-(CM(2,5)) (CM_IVD(1,3))-(CM(2,6))];
 %superior point on the vert and inferior on the ivd
 L_no_r(2).VT = [L_no_r(2).VT(:,1)+t(1) L_no_r(2).VT(:,2)+t(2) L_no_r(2).VT(:,3)+t(3)];

 CM(2,1:3) = [CM(2,1)+t(1) CM(2,2)+t(2) CM(2,3)+t(3)];
 CM(2,4:6) = [CM(2,4)+t(1) CM(2,5)+t(2) CM(2,6)+t(3)];
 CM(2,7:9) = [CM(2,7)+t(1) CM(2,8)+t(2) CM(2,9)+t(3)];

% plot3(L_no_r(2).VT(:,1),L_no_r(2).VT(:,2),L_no_r(2).VT(:,3),'.y'),hold on;

% -------- vert2 -IVD2
t = [ (CM(2,1))-(CM_IVD(2,4)) (CM(2,2))-(CM_IVD(2,5)) (CM(2,3))-(CM_IVD(2,6))]; %superior point on the vert and inferior on the ivd
         
mesh_struct_IVD2(2).V =  [mesh_struct_IVD2(2).V(:,1)+(t(1)) ...
  mesh_struct_IVD2(2).V(:,2)+(t(2)) mesh_struct_IVD2(2).V(:,3)+(t(3)) ];
 
CM_IVD(2,1:3) = [CM_IVD(2,1)+(t(1)) CM_IVD(2,2)+(t(2)) CM_IVD(2,3)+(t(3))];
CM_IVD(2,4:6) = [CM_IVD(2,4)+(t(1)) CM_IVD(2,5)+(t(2)) CM_IVD(2,6)+(t(3))];
CM_IVD(2,7:9) = [CM_IVD(2,7)+(t(1)) CM_IVD(2,8)+(t(2)) CM_IVD(2,9)+(t(3))];


% plot3(mesh_struct_IVD2(2).V(:,1),mesh_struct_IVD2(2).V(:,2),mesh_struct_IVD2(2).V(:,3),'.b'),hold on;

% -------- IVD2 - vert3

 t = [ (CM_IVD(2,1))-(CM(3,4)) (CM_IVD(2,2))-(CM(3,5)) (CM_IVD(2,3))-(CM(3,6))];
 %superior point on the vert and inferior on the ivd
 L_no_r(3).VT = [L_no_r(3).VT(:,1)+t(1) L_no_r(3).VT(:,2)+t(2) L_no_r(3).VT(:,3)+t(3)];

 CM(3,1:3) = [CM(3,1)+t(1) CM(3,2)+t(2) CM(3,3)+t(3)];
 CM(3,4:6) = [CM(3,4)+t(1) CM(3,5)+t(2) CM(3,6)+t(3)];
 CM(3,7:9) = [CM(3,7)+t(1) CM(3,8)+t(2) CM(3,9)+t(3)];

% plot3(L_no_r(3).VT(:,1),L_no_r(3).VT(:,2),L_no_r(3).VT(:,3),'.y'),hold on;
 % -------- vert3 -IVD3
t = [ (CM(3,1))-(CM_IVD(3,4)) (CM(3,2))-(CM_IVD(3,5)) (CM(3,3))-(CM_IVD(3,6))]; %superior point on the vert and inferior on the ivd
         
mesh_struct_IVD2(3).V =  [mesh_struct_IVD2(3).V(:,1)+(t(1)) ...
  mesh_struct_IVD2(3).V(:,2)+(t(2)) mesh_struct_IVD2(3).V(:,3)+(t(3)) ];
 
CM_IVD(3,1:3) = [CM_IVD(3,1)+(t(1)) CM_IVD(3,2)+(t(2)) CM_IVD(3,3)+(t(3))];
CM_IVD(3,4:6) = [CM_IVD(3,4)+(t(1)) CM_IVD(3,5)+(t(2)) CM_IVD(3,6)+(t(3))];
CM_IVD(3,7:9) = [CM_IVD(3,7)+(t(1)) CM_IVD(3,8)+(t(2)) CM_IVD(3,9)+(t(3))];

%plot3(mesh_struct_IVD2(3).V(:,1),mesh_struct_IVD2(3).V(:,2),mesh_struct_IVD2(3).V(:,3),'.b'),hold on;

% -------- IVD3 - vert4

 t = [ (CM_IVD(3,1))-(CM(4,4)) (CM_IVD(3,2))-(CM(4,5)) (CM_IVD(3,3))-(CM(4,6))];
 %superior point on the vert and inferior on the ivd
 L_no_r(4).VT = [L_no_r(4).VT(:,1)+t(1) L_no_r(4).VT(:,2)+t(2) L_no_r(4).VT(:,3)+t(3)];

 CM(4,1:3) = [CM(4,1)+t(1) CM(4,2)+t(2) CM(4,3)+t(3)];
 CM(4,4:6) = [CM(4,4)+t(1) CM(4,5)+t(2) CM(4,6)+t(3)];
 CM(4,7:9) = [CM(4,7)+t(1) CM(4,8)+t(2) CM(4,9)+t(3)];

% plot3(L_no_r(4).VT(:,1),L_no_r(4).VT(:,2),L_no_r(4).VT(:,3),'.y'),hold on;
 
  % -------- vert4 -IVD4
t = [ (CM(4,1))-(CM_IVD(4,4)) (CM(4,2))-(CM_IVD(4,5)) (CM(4,3))-(CM_IVD(4,6))]; %superior point on the vert and inferior on the ivd
         
mesh_struct_IVD2(4).V =  [mesh_struct_IVD2(4).V(:,1)+(t(1)) ...
  mesh_struct_IVD2(4).V(:,2)+(t(2)) mesh_struct_IVD2(4).V(:,3)+(t(3)) ];
 
CM_IVD(4,1:3) = [CM_IVD(4,1)+(t(1)) CM_IVD(4,2)+(t(2)) CM_IVD(4,3)+(t(3))];
CM_IVD(4,4:6) = [CM_IVD(4,4)+(t(1)) CM_IVD(4,5)+(t(2)) CM_IVD(4,6)+(t(3))];
CM_IVD(4,7:9) = [CM_IVD(4,7)+(t(1)) CM_IVD(4,8)+(t(2)) CM_IVD(4,9)+(t(3))];

%plot3(mesh_struct_IVD2(4).V(:,1),mesh_struct_IVD2(4).V(:,2),mesh_struct_IVD2(4).V(:,3),'.b'),hold on;

% -------- IVD4 - vert5

 t = [ (CM_IVD(4,1))-(CM(5,4)) (CM_IVD(4,2))-(CM(5,5)) (CM_IVD(4,3))-(CM(5,6))];
 %superior point on the vert and inferior on the ivd
 L_no_r(5).VT = [L_no_r(5).VT(:,1)+t(1) L_no_r(5).VT(:,2)+t(2) L_no_r(5).VT(:,3)+t(3)];

 CM(5,1:3) = [CM(5,1)+t(1) CM(5,2)+t(2) CM(5,3)+t(3)];
 CM(5,4:6) = [CM(5,4)+t(1) CM(5,5)+t(2) CM(5,6)+t(3)];
 CM(5,7:9) = [CM(5,7)+t(1) CM(5,8)+t(2) CM(5,9)+t(3)];

% plot3(L_no_r(5).VT(:,1),L_no_r(5).VT(:,2),L_no_r(5).VT(:,3),'.y'),hold on;

%{
for(j=1:2)
plot3(CM(j,1),CM(j,2),CM(j,3),'Og'),hold on; %down
plot3(CM(j,4),CM(j,5),CM(j,6),'Og'),hold on; %up
end

for(k=1:1)
plot3(CM_IVD(k,1),CM_IVD(k,2),CM_IVD(k,3),'og'),hold on; %bottom
plot3(CM_IVD(k,4),CM_IVD(k,5),CM_IVD(k,6),'og'),hold on; %up
end
%}
end