%{
New script to give the orientation to the Lumbar spine

I need:
- distance between the sup and inf surface of each body
- boundary conditions to evaluate the min distance

%}

function  [CM, CM_IVD,mesh_struct_IVD2,Lmes]  = orientation_gui( Lmes, mesh_struct_IVD2, nbodies, alpha, IVD, heighV )
%name = sprintf('./Lmes_v2.mat');
%load(name);

    
    %plot3(Lmes(1).VT(:,1),Lmes(1).VT(:,2),Lmes(1).VT(:,3),'.b'),hold on;
    %plot3(Lmes(2).VT(:,1),Lmes(2).VT(:,2),Lmes(2).VT(:,3),'.b'),hold on;
    %plot3(Lmes(3).VT(:,1),Lmes(3).VT(:,2),Lmes(3).VT(:,3),'.g'),hold on;
   % plot3(mesh_struct_IVD2(1).V(:,1),mesh_struct_IVD2(1).V(:,2),mesh_struct_IVD2(1).V(:,3),'.b'),hold on




alpha_lum = alpha;
alpha_lum_r = pi*alpha_lum/180; %angle in radians
hl=0; %   total length
%for the vertebrae
    for(j =1:5)
        %plot3( Lmes(j).VT(:,1),Lmes(j).VT(:,2),Lmes(j).VT(:,3),'.b'), hold on;


        logicFace=Lmes(j).Cb==2; % select only the face outside
        Face_ext= Lmes(j).Fb(logicFace,:);

        [N]=patchNormal(Face_ext, Lmes(j).VT); 
        faceBoundaryMarker=zeros(size(Face_ext,1),1);

        faceBoundaryMarker(N(:,1)<-0.5)=1; %Left
        faceBoundaryMarker(N(:,1)>0.5)=2; %Right
        faceBoundaryMarker(N(:,2)<-0.5)=3; %Front
        faceBoundaryMarker(N(:,2)>0.5)=4; %Back
        faceBoundaryMarker(N(:,3)<-0.5)=5; %Top
        faceBoundaryMarker(N(:,3)>0.5)=6; %Bottom

        logicFace2=faceBoundaryMarker==5;
        F(j).vert_sup= Face_ext(logicFace2,:);
        bc(j).vert_sup=unique(F(j).vert_sup(:));
        % find nodes associated to the faces
        bc(j).V_sup = Lmes(j).VT(bc(j).vert_sup,:);
        bc(j).V_bod_sup = bc(j).V_sup(bc(j).V_sup(:,2)>-20,:); %Superior surface only v body 

        %x & y centre 
        mean_xy = [ mean(bc(j).V_bod_sup(:,1)) mean(bc(j).V_bod_sup(:,2))];
        %find the mean z
        centre = bc(j).V_bod_sup(find(bc(j).V_bod_sup(:,1)<mean_xy(1)+1.5 & bc(j).V_bod_sup(:,1)>mean_xy(1)-1.5 &...
            bc(j).V_bod_sup(:,2)<mean_xy(2)+1.5  & bc(j).V_bod_sup(:,2)>mean_xy(2)-1.5),:) ;
        %plot3(centre(1),centre(2),centre(3),'Or'),hold on;
        % centroid of the endplate    
        bc(j).V_cent_sup =[mean_xy mean(centre(:,3))];


        logicFace2=faceBoundaryMarker==6;
        F(j).vert_inf=Face_ext(logicFace2,:);
        bc(j).vert_inf=unique(F(j).vert_inf(:));
        % find nodes associated to the faces
        bc(j).V_inf = Lmes(j).VT(bc(j).vert_inf,:);
        bc(j).V_bod_inf = bc(j).V_inf(bc(j).V_inf(:,2)>-20,:);%inferior surface only v body 

        %x & y centre 
        mean_xy = [ mean(bc(j).V_bod_inf(:,1)) mean(bc(j).V_bod_inf(:,2))];
        %find the mean z
        centre = bc(j).V_bod_inf(find(bc(j).V_bod_inf(:,1)<mean_xy(1)+1.5 & bc(j).V_bod_inf(:,1)>mean_xy(1)-1.5 &...
        bc(j).V_bod_inf(:,2)<mean_xy(2)+1.5  & bc(j).V_bod_inf(:,2)>mean_xy(2)-1.5),:) ;
        % centroid of the endplate    
        bc(j).V_cent_inf =[mean_xy mean(centre(:,3))];
        
        z = [ bc(j).V_cent_sup(:,3) bc(j).V_cent_inf(:,3) ] ; %z centre
        CG(j,:) = [mean_xy, mean(z)];
        %plot3( bc(j).V_bod_inf(:,1),bc(j).V_bod_inf(:,2),bc(j).V_bod_inf(:,3),'.r'), hold on;
        %plot3( bc(j).V_bod_sup(:,1),bc(j).V_bod_sup(:,2),bc(j).V_bod_sup(:,3),'.r'),hold on;
        %plot3(bc(j).V_cent_inf(:,1),bc(j).V_cent_inf(:,2),bc(j).V_cent_inf(:,3),'Ob'),hold on;
        %plot3(bc(j).V_cent_sup(:,1),bc(j).V_cent_sup(:,2),bc(j).V_cent_sup(:,3),'Ob'),hold on;

        % 2.  total length
        hl = hl + (bc(j).V_cent_sup(:,3)-bc(j).V_cent_inf(:,3));


    end

    
    %for the IVD (mesh_struct_IVD2)
    for(j =1:4)
        
        % 1 top 2 bottom
       % boundary facets meshStruct_IVD2.Fb=Fb;
        
       %--- superior surface
        logicFace=mesh_struct_IVD2(j).faceBoundaryMarker==2; % select only the face outside
        I(j).vert_sup= mesh_struct_IVD2(j).Fb(logicFace,:);
        
        bcIVD(j).vert_sup=unique(I(j).vert_sup(:));
        % find nodes associated to the faces
        bcIVD(j).V_sup =mesh_struct_IVD2(j).V(bcIVD(j).vert_sup,:); 
       % plot3(mesh_struct_IVD2(j).V(:,1),mesh_struct_IVD2(j).V(:,2),mesh_struct_IVD2(j).V(:,3),'.r'),hold on;
        %plot3(bcIVD(j).V_sup(:,1),bcIVD(j).V_sup(:,2),bcIVD(j).V_sup(:,3),'Oy'),hold on;
        
        %--- inferior surface
        logicFace2=mesh_struct_IVD2(j).faceBoundaryMarker==1; % select only the face outside
        I(j).vert_inf= mesh_struct_IVD2(j).Fb(logicFace2,:);
        bcIVD(j).vert_inf =unique(I(j).vert_inf(:));
        % find nodes associated to the faces
        bcIVD(j).V_inf =mesh_struct_IVD2(j).V(bcIVD(j).vert_inf,:);
        
        %plot3(mesh_struct_IVD2(j).V(:,1),mesh_struct_IVD2(j).V(:,2),mesh_struct_IVD2(j).V(:,3),'.g'),hold on;
       % plot3(bcIVD(j).V_inf(:,1),bcIVD(j).V_inf(:,2),bcIVD(j).V_inf(:,3),'Og'),hold on;
        
        
              
        % centre sup surface
        centreIVD_sup = mean(bcIVD(j).V_sup);
        bcIVD(j).V_centre_sup = centreIVD_sup;
        %plot3(centreIVD_sup(1),centreIVD_sup(2),centreIVD_sup(3),'Ob'),hold on;

        % centre inf surface
        centreIVD_inf = mean(bcIVD(j).V_inf);
         bcIVD(j).V_centre_inf = centreIVD_inf;
        %plot3(centreIVD_inf(1),centreIVD_inf(2),centreIVD_inf(3),'Or'),hold on;
          
        % centre Disc
        CG_IVD(j,:) = mean(mesh_struct_IVD2(j).V);
        %plot3( CG_IVD(:,1),CG_IVD(:,2),CG_IVD(:,3),'Or'), hold on;
        
    end
    
c = hl;% + sum(IVD); %(mm)
%radius corresponding circumference
r = c /( 2*sin(alpha_lum_r) );

y_cent = 0;
z_cent = 0;
th = 0: alpha_lum_r/100: alpha_lum_r;
yunit = r * cos(th) + y_cent;
zunit = r * sin(th) + z_cent;
x = zeros(length(zunit)); 
%figure;
%plot3(x,yunit, zunit, '-.g'),hold on;


%-------  traslate
% starting from the 1st vertyebrae, I can move the upper point of the
% vertebra to be coincident with the first point on the circumference
alp = 0;
y_up = yunit(1);
z_up = zunit(1);
eps = 0;
vert_r = [];
count_v = 1;
count_i = 1;
count_colv = 1;
count_coli = 1;
count_ivd = 1;
z_rec = 0; y_rec = 0;
CM_curve_y = []; CM_curve_z = [];
ii = 0;
    for(body = 1:9)
        %body;

       if(mod(body,2) ~= 0) 
            % vertebrae
                    %moving only on y and z
                 
    %plot3(Lmes(count_v).VT(:,1),Lmes(count_v).VT(:,2),Lmes(count_v).VT(:,3),'.g'),hold on;
             t = [ bc(count_v).V_cent_sup(2)-y_up, bc(count_v).V_cent_sup(3)-z_up ];

    % move the dataset
            %traslation

            Lmes(count_v).VT = [ Lmes(count_v).VT(:,1), Lmes(count_v).VT(:,2)-t(1), Lmes(count_v).VT(:,3)-abs(t(2)) ];


            % update the CM vector
            bc(count_v).V_cent_sup(:) = [ bc(count_v).V_cent_sup(1), bc(count_v).V_cent_sup(2)-(t(1)), bc(count_v).V_cent_sup(3)-abs(t(2)) ];
            bc(count_v).V_cent_inf(:) = [ bc(count_v).V_cent_inf(1), bc(count_v).V_cent_inf(2)-(t(1)), bc(count_v).V_cent_inf(3)-abs(t(2)) ];
            CG(count_v,:) = [CG(count_v,1),CG(count_v,2)-t(1),CG(count_v,3)-abs(t(2))];

            % evaluate the vector between the points on the upper and bottom surfaces
            %v1 vector (I m working on y and z) defined between the vertebrae
            v1 = [  bc(count_v).V_cent_inf(2)- bc(count_v).V_cent_sup(2),  bc(count_v).V_cent_inf(3)- bc(count_v).V_cent_sup(3)];
            %find the corrispondent point on the circumference


            d = sqrt( (bc(count_v).V_cent_sup(2)-bc(count_v).V_cent_inf(2))^2+ (bc(count_v).V_cent_sup(3)-bc(count_v).V_cent_inf(3))^2 ); %corda
          
            alp = asin(d/(2*abs(r)))*2+alp;
            %plot(y_up, z_up,'*m');
            y_rec = r * cos(alp) + y_cent; % +y_rec;
            z_rec = r * sin(alp) + z_cent;% +z_rec;
            %plot(y_rec, z_rec, '^B'),hold on;

            %vector between the points belonging to the circumference

            v2 = [ y_rec-y_up, z_rec-z_up ];

            %find the angle in between them and rotate the body of that amount
            dot_prod = dot(v1,v2);
            angle = -acos(dot_prod/(norm(v1)*norm(v2)));
            mat = [cos(angle),-sin(angle);sin(angle),cos(angle)];
            %////////////////////////// ok
           CM(count_v,1:6) = [bc(count_v).V_cent_inf, bc(count_v).V_cent_sup];
            CM(count_v,7:9) = CG(count_v,:);

           %size(CM)
            [ Lmes(count_v).VT , CM(count_v,:)] = ...
                moving_rot_post_mesh_gui(Lmes(count_v).VT, CM(count_v,:), t, mat);
           
            rot_angle(count_v) = angle; 
            
            count_v = count_v+1 ;
       
       else
            
      
            %////////////////////////// ok
           CM_IVD(count_ivd,1:6) = [bcIVD(count_ivd).V_centre_inf, bcIVD(count_ivd).V_centre_sup];
           CM_IVD(count_ivd,7:9) = CG_IVD(count_ivd,:);
           
            %plot3(CM_IVD(count_ivd,1),CM_IVD(count_ivd,2),CM_IVD(count_ivd,3),'Og'),hold on;
            %plot3(mesh_struct_IVD2(count_ivd).V(:,1),mesh_struct_IVD2(count_ivd).V(:,2),mesh_struct_IVD2(count_ivd).V(:,3),'.m' ),hold on;
            [ mesh_struct_IVD2(count_ivd).V , CM_IVD(count_ivd,:)] = ...
                moving_rot_post_mesh_gui(mesh_struct_IVD2(count_ivd).V, CM_IVD(count_ivd,:), t, mat);
            % plot3(CM_IVD(count_ivd,1),CM_IVD(count_ivd,2),CM_IVD(count_ivd,3),'Og'),hold on;
            %plot3(mesh_struct_IVD2(count_ivd).V(:,1),mesh_struct_IVD2(count_ivd).V(:,2),mesh_struct_IVD2(count_ivd).V(:,3),'.m' ),hold on;
        
            
            count_ivd = count_ivd +1;
           
       end
       
    end
      


 [CM, CM_IVD,mesh_struct_IVD2,Lmes] = position_curve_gui2(Lmes, mesh_struct_IVD2, CM, CM_IVD, alpha, hl,rot_angle, IVD, heighV);  
  
   
end