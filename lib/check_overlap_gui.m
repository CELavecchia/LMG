%function to remove the overlap. L_body is the main body and L_ped_dx is
%the body where the points have to be removed

function ped_new =check_overlap( L_body, L_ped_dx )


%remove overlapping points on the pedicles, and on the lamina

%first at all, initial cleaning
L_bod_sel = L_body(L_body(:,3)>mean(L_body(:,3)),:);
L_ped_dx = L_ped_dx( L_ped_dx(:,2)<= min(L_bod_sel(:,2))+2,:);



%check points and check the z values
eps = 0.3;
count = 1; count_new = 1 ;
[r,c] = find(L_body(:,2)>=min(L_ped_dx(:,2)));
if(isempty(r)==0)
selection_L_body= L_body(r,:);
clear r; clear c;
[r,c] = find(selection_L_body(:,3)>=min(L_ped_dx(:,3)));
selection_L_body2 = selection_L_body(r,:);
%plot3(selection_L_body2(:,1),selection_L_body2(:,2),selection_L_body2(:,3),'.b'),hold on;

ped_new = []; ped_save = [];
for(j = 1:size(selection_L_body2,1))
     for(k = 1:size(L_ped_dx,1))
        if( L_ped_dx(k,2)< selection_L_body2(j,2)+eps && L_ped_dx(k,2)>= selection_L_body2(j,2)-eps && ...
            L_ped_dx(k,1)< selection_L_body2(j,1)+eps && L_ped_dx(k,1)>= selection_L_body2(j,1)-eps && ...
            L_ped_dx(k,3)< selection_L_body2(j,3)+eps && L_ped_dx(k,3)>= selection_L_body2(j,3)-eps)
         %clean better
        ped_save(count,:) = L_ped_dx(k,:) ; 
        count = count+1;
        
        end 
     end   
end
%plot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'.r'),hold on;
%plot3(ped_save(:,1),ped_save(:,2),ped_save(:,3),'*g'),hold on;  %check it

ped_extra = L_ped_dx(L_ped_dx(:,2)<mean(ped_save(:,2)),:);
ped_save = [ ped_save;ped_extra];

%plot3(ped_save(:,1),ped_save(:,2),ped_save(:,3),'Oy'),hold on;  %check it
%plot3(L_ped_dx(:,1),L_ped_dx(:,2),L_ped_dx(:,3),'Oy'),hold on;  %check it

%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.b'),hold on;
ped_new = ped_save;
%ped_new = setdiff(L_ped_dx,ped_save,'rows');


%plot3(ped_new(:,1), ped_new(:,2), ped_new(:,3),'*r'),hold on;
end

end