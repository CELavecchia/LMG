%check the alignment between the vertebrae and IVD

function [CM, CM_IVD, mesh_struct_IVD2, L_no_r ] =check_position_arch_gui(CM, CM_IVD, mesh_struct_IVD2, L_no_r)


%align the inferior point of the IVD with the superior of the following
%vertebrae and the inferior point of the vertebrae with the superior ones
%of the IVD
%figure(3)
count_ivd = 1;
count_v =1;
%{
for(j=1:5)
plot3(L_no_r(j).VT(:,1),L_no_r(j).VT(:,2),L_no_r(j).VT(:,3),'.r'),hold on;
end
for(k=1:4)
         plot3(mesh_struct_IVD2(k).V(:,1),mesh_struct_IVD2(k).V(:,2),mesh_struct_IVD2(k).V(:,3),'.r'),hold on;
end
%}
figure
plot3(L_no_r(1).VT(:,1),L_no_r(1).VT(:,2),L_no_r(1).VT(:,3),'.r'),hold on;

for (body=1:8)
    
     if(mod(body,2) ~= 0) %if vert body
        
         t = [ (CM(count_v,1))-(CM_IVD(count_ivd,4)) (CM(count_v,2))-(CM_IVD(count_ivd,5)) (CM(count_v,3))-(CM_IVD(count_ivd,6))] %inferior point on the vert and superior on the ivd
         
         mesh_struct_IVD2(count_ivd).V =  [mesh_struct_IVD2(count_ivd).V(:,1)+t(1) ...
             mesh_struct_IVD2(count_ivd).V(:,2)-t(2) mesh_struct_IVD2(count_ivd).V(:,3)+t(3) ];
         CM_IVD(count_ivd,1:3) = [CM_IVD(count_ivd,1)+t(1) CM_IVD(count_ivd,2)-t(2) CM_IVD(count_ivd,3)+t(3)];
         CM_IVD(count_ivd,4:6) = [CM_IVD(count_ivd,4)+t(1) CM_IVD(count_ivd,5)-t(2) CM_IVD(count_ivd,6)+t(3)];
         CM_IVD(count_ivd,7:9) = [CM_IVD(count_ivd,7)+t(1) CM_IVD(count_ivd,8)-t(2) CM_IVD(count_ivd,9)+t(3)];

        
         plot3(mesh_struct_IVD2(count_ivd).V(:,1),mesh_struct_IVD2(count_ivd).V(:,2),mesh_struct_IVD2(count_ivd).V(:,3),'.r'),hold on;
         count_v = count_v+1;
    
     else %if disc
        
          
         t = [ abs(CM_IVD(count_ivd,1))-abs(CM(count_v,4)) abs(CM_IVD(count_ivd,2))-abs(CM(count_v,5)) abs(CM_IVD(count_ivd,3))-abs(CM(count_v,6))]
         %superior point on the vert and inferior on the ivd
         L_no_r(count_v).VT = [L_no_r(count_v).VT(:,1)+t(1) L_no_r(count_v).VT(:,2)+t(2) L_no_r(count_v).VT(:,3)+t(3)];
        
         CM(count_v,1:3) = [CM(count_v,1)+t(1) CM(count_v,2)-t(2) CM(count_v,3)+t(3)];
         CM(count_v,4:6) = [CM(count_v,4)+t(1) CM(count_v,5)-t(2) CM(count_v,6)+t(3)];
         CM(count_v,7:9) = [CM(count_v,7)+t(1) CM(count_v,8)-t(2) CM(count_v,9)+t(3)];
         
         plot3(L_no_r(count_v).VT(:,1),L_no_r(count_v).VT(:,2),L_no_r(count_v).VT(:,3),'.y'),hold on;


        count_ivd = count_ivd +1;
     end
         
    

end
%}

%figure(3)
for(j=1:5)
plot3(CM(j,1),CM(j,2),CM(j,3),'or'),hold on; %down
plot3(CM(j,4),CM(j,5),CM(j,6),'or'),hold on; %up
end

for(k=1:4)
plot3(CM_IVD(k,1),CM_IVD(k,2),CM_IVD(k,3),'og'),hold on; %bottom
plot3(CM_IVD(k,4),CM_IVD(k,5),CM_IVD(k,6),'og'),hold on; %up
end

end