%{ 
points on the ellipse
INPUT : - center coordinates - semiaxes
OUTPUT : - P POINTS on the ellipse
         - p_x, p_z gridded points

add the centre and line orientation

%}


function [ L_ped ] = ellipse_gui( P_sag, P_trav, sc_d, sc_w, b, a, angl_sagit, angl_trasv )

    %first line, sagittal plane [Xs, Ys] 
    ms = atan(pi/2-angl_sagit); qs = P_sag(2)-ms*P_sag(1);
    x = (P_sag(1) - sc_w/2):0.5:0;
    
    y = ms*x+qs;
    %second line, trasversal plane
    mt = atan(angl_trasv); qt = P_trav(3)-mt*P_trav(2);
    z = mt*y+qt;
   % plot3(x,y,z,'.b'),hold on;
    
    

    t=-pi:0.1:pi;
    count =1;
    for(j = 1:length(y) )
    
        x1(count:count+length(t)-1)=x(j)+a/2*cos(t);
        z1(count:count+length(t)-1)=z(j)+b/2*sin(t);
        y1(count:count+length(t)-1)= ones(size(t,1))*y(j);
        count = count + length(t);
        
        
    end
       
          L_ped = [ x1', y1', z1'; -x1',y1',z1'];
          %plot3(L_ped(:,1),L_ped(:,2),L_ped(:,3),'.r'),hold on;
          


end
