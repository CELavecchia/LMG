function [ mesh_struct_IVD] = enforce_contacts(mesh_struct_IVD,EP)

    warning off
    
    %import 
    
    %pre-processing vertebral upper and lower curves 
    interpFunction = fitting_surfaces_gui(EP);
    
    
    n_layer = 9;
    n_points = 9;
    count = 1;
   
    %V2 =reshape(mesh_struct_IVD(1).V,
    %[size(mesh_struct_IVD(1).V(:,:),1)/n_layer, 3, n_layer]); % I don t
    %understand why it s not woeking
   
    
    %fitting the upper and lower surfaces and interpolate the nodes
    for(w=1:4)
    
        for j = 1:(size(mesh_struct_IVD(w).V,1)/n_layer) 
            ind =j;        
                       
           for(k =1:n_layer)
                
                mesh_struct_line(k,:) = mesh_struct_IVD(w).V(ind,:);
                %save all the elements in the same line
                ind_vect(k) = ind;
                 ind = ind + size(mesh_struct_IVD(w).V,1)/n_layer;
            end
            
            point_up = interpFunction(w).inf( mesh_struct_line(1,1), mesh_struct_line(1,2));
            point_bot = interpFunction(w+1).sup( mesh_struct_line(end,1), mesh_struct_line(end,2));
            P = [mesh_struct_line(1,1), mesh_struct_line(1,2) point_up; mesh_struct_line(end,1), mesh_struct_line(end,2) point_bot];
            
            mesh_struct_line = evenlySampleCurve(P, n_points,'pchip',0);
             
            mesh_struct_IVD(w).V(ind_vect,:) = mesh_struct_line;
        
            
        end
    end
   % V3 = reshape(V3, size(V3,1)*n_layer,3,1);
%    figure
%    plot3(mesh_struct_IVD(1).V(:,1),mesh_struct_IVD(1).V(:,2),mesh_struct_IVD(1).V(:,3),'.r'),hold on
%    plot3(mesh_struct_line(:,1),mesh_struct_line(:,2),mesh_struct_line(:,3),'ob')
%    plot3(mesh_struct_IVD(1).V(:,1),mesh_struct_IVD(1).V(:,2),mesh_struct_IVD(1).V(:,3),'.g'),hold on
%     

end