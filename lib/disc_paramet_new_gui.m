function   disc_nodes = disc_paramet_new_gui( vertices, EPWu_half, EPDu, IVD);


 
    %evaluate the dimensions and scaling factor
    
    % width
    max_x = max(vertices(:,1));
    min_x = min(vertices(:,1));
    width_h = (max_x - min_x)/2; %half width
    sc_aw = EPWu_half/width_h; %scaling factor
   % plot3(anulus_nodes(:,1)*sc_aw,anulus_nodes(:,2),anulus_nodes(:,3),'.r'),hold on, grid on;
    %plot3(nucleus_nodes(:,1)*sc_aw,nucleus_nodes(:,2),nucleus_nodes(:,3),'.r'),hold on, grid on;
    vertices_w = [ vertices(:,1)*sc_aw, vertices(:,2), vertices(:,3)];
    
    %depth
    anulus_meanx = vertices_w(vertices_w(:,1)<1 & vertices_w(:,1)>-1,:);
   % plot3(anulus_meanx(:,1),anulus_meanx(:,2),anulus_meanx(:,3),'.g'),hold on;
    max_y = max(anulus_meanx(:,2));
    min_y = min(anulus_meanx(:,2));
    depth = (max_y - min_y); %half depth
    sc_d = EPDu/depth; %scaling factor
    %plot3(anulus_nodes(:,1),anulus_nodes(:,2)*sc_d,anulus_nodes(:,3),'.r'),hold on, grid on;
    %plot3(nucleus_nodes(:,1),nucleus_nodes(:,2)*sc_d,nucleus_nodes(:,3),'.r'),hold on, grid on;
    disc_nodes = [vertices_w(:,1),vertices_w(:,2)*sc_d,vertices_w(:,3)];
   
    %height
     nucleus_meanz = disc_nodes( disc_nodes(:,1)<1 & disc_nodes(:,1)>-1 ...
         & disc_nodes(:,2)<5 & disc_nodes(:,2)>-5,:);
   % plot3(nucleus_nodes(:,1),nucleus_nodes(:,2),nucleus_nodes(:,3),'.g'),hold on;
    max_z = max(nucleus_meanz(:,3));
    min_z = min(nucleus_meanz(:,3));
    height = (max_z - min_z); %height in the central part of the disc
    sc_z = IVD/height; %scaling factor
    
   disc_nodes = [disc_nodes(:,1),disc_nodes(:,2),disc_nodes(:,3)*sc_z];
   
    %plot3(anulus_nodes(:,1),anulus_nodes(:,2),anulus_nodes(:,3),'.r'),hold on, grid on;
    %plot3(nucleus_nodes(:,1),nucleus_nodes(:,2),nucleus_nodes(:,3),'.r'),hold on, grid on;


     
end