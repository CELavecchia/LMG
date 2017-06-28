%{
New script to give the orientation to the Lumbar spine

I need:
- distance between the sup and inf surface of each body
- boundary conditions to evaluate the min distance

%}
clc; clear all; close all;

function orientation_post_mesh_gui( Lmes )
%name = sprintf('./Lmes_v2.mat');
%load(name);





alpha_lum = 43.49;
alpha_lum_r = pi*alpha_lum/180; %angle in radians
hl=0; %   total length

for(j =1:5)
    plot3( Lmes(j).VT(:,1),Lmes(j).VT(:,2),Lmes(j).VT(:,3),'.b'), hold on;
    
    
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
    
    %plot3( bc(j).V_bod_inf(:,1),bc(j).V_bod_inf(:,2),bc(j).V_bod_inf(:,3),'.r'), hold on;
    %plot3( bc(j).V_bod_sup(:,1),bc(j).V_bod_sup(:,2),bc(j).V_bod_sup(:,3),'.r'),hold on;
    %plot3(bc(j).V_cent_inf(:,1),bc(j).V_cent_inf(:,2),bc(j).V_cent_inf(:,3),'Ob'),hold on;
    %plot3(bc(j).V_cent_sup(:,1),bc(j).V_cent_sup(:,2),bc(j).V_cent_sup(:,3),'Ob'),hold on;
   
    % 2.  total length
    hl = hl + (bc(j).V_cent_sup(:,3)-bc(j).V_cent_inf(:,3));
    

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
for(body = 1:5)
    %body;
    
  %  if(mod(body,2) ~= 0) 
        % vertebrae
                %moving only on y and z
        %ii = ii+1;
         t = [ bc(body).V_cent_sup(2)-y_up, bc(body).V_cent_sup(3)-z_up ]
        
% move the dataset
        %traslation
     
        Lmes(body).VT = [ Lmes(body).VT(:,1), Lmes(body).VT(:,2)-t(1), Lmes(body).VT(:,3)-abs(t(2)) ];
        
        
        % update the CM vector
        bc(body).V_cent_supf(:) = [ bc(body).V_cent_sup(1), bc(body).V_cent_sup(2)-abs(t(1)), bc(body).V_cent_sup(3)-abs(t(2)) ];
        
           
        % evaluate the vector between the points on the upper and bottom surfaces

        %v1 vector (I m working on y and z) defined between the vertebrae
        v1 = [  bc(body).V_cent_inf(2)- bc(body).V_cent_sup(2),  bc(body).V_cent_inf(3)- bc(body).V_cent_sup(3)];
        %find the corrispondent point on the circumference
        

        d = sqrt( (bc(body).V_cent_sup(2)-bc(body).V_cent_sup(2))^2+ (bc(body).V_cent_supf(3)-bc(body).V_cent_sup(3))^2 ); %corda
       %{
        alp_test = asin( ( d+eps )/(2*abs(r)))*2+alp;
        y_rec2 = r * cos(alp_test) ;
        z_rec2 = r * sin(alp_test) ;
        plot3(0,y_rec2, z_rec2, 'Og'),hold on; % ending of the body
      %}

        alp = asin(d/(2*abs(r)))*2+alp;
        
        %plot(y_up, z_up,'*m');
        
        y_rec = r * cos(alp) + y_cent; % +y_rec;
        z_rec = r * sin(alp) + z_cent;% +z_rec;
        %plot(y_rec, z_rec, '^B'),hold on;
        
        %vector between the points belonging to the circumference
     
        v2 = [ y_rec-y_up, z_rec-z_up ];
       
        %find the angle in between them and rotate the body of that amount
        dot_prod = dot(v1,v2);
        angle = acos(dot_prod/(norm(v1)*norm(v2)));
        mat = [cos(angle),-sin(angle);sin(angle),cos(angle)];
        %////////////////////////// ok
       CM = [bc(body).V_cent_inf, bc(body).V_cent_sup];
        
        [ Lmes(body).VT , CM] = ...
            moving_rot_post_mesh(Lmes(body).VT, CM, count_colv, count_v, t, mat);
               
        plot3(Lmes(body).VT(:,1),Lmes(body).VT(:,2),Lmes(body).VT(:,3),'.r'),hold on;
       
        %plot(CM(1,5), CM(1,6), 'Oy');
        inf = [ CM(1), CM(count_v,2), CM(count_v,3) ];
       % plot3(inf(1),inf(2),inf(3),'Om');
        count_v = count_v+1 ;
        count_colv = count_colv+3 ;
        
end
end