%function to remove the overlap. L_body is the main body and L_ped_dx is
%the body where the points have to be removed

function ped_new =check_overlap_ped_pro( L_ped, L_proc_dx )


%remove overlapping points on the pedicles, and on the lamina

%check it
L_proc_selec = L_proc_dx(L_proc_dx(:,2)>=min(L_ped(:,2)),:);
L_ped_selec = L_ped(L_ped(:,2)<=max(L_proc_dx(:,2)),:);

%plot3(L_proc_selec(:,1),L_proc_selec(:,2),L_proc_selec(:,3),'.r'),hold on;
%plot3(L_ped_selec(:,1),L_ped_selec(:,2),L_ped_selec(:,3),'Ob');

%clean it a little bit better
eps=0.05; count =1;
size(L_ped_selec);
size(L_proc_selec);

%between the L_ped_selec, evaluate the distance from the
%min(L_proc_selec(:,1) and save only those that are further than 0.5
[max_proc,ind] = max(L_proc_selec(:,2));
max_proc = L_proc_selec(ind,:);
%{
eps = 0.4;
count =1;
for(j = 1: size(L_ped_selec))
   for(k = 1:size(L_proc_selec))
    if( L_ped_selec(j,1)<L_proc_selec(k,1)+eps & ...
            L_ped_selec(j,1)>L_proc_selec(k,1)-eps)
        remove(count,:) = L_ped_selec(j,:);
        count = count+1;
    end
   end
end
%}

ped_new = setdiff(L_ped,L_ped_selec,'rows');
%plot3(ped_new(:,1),ped_new(:,2),ped_new(:,3),'.g'),hold on;
%{
L_ped_save = setdiff(L_ped,remove,'rows');
plot3(L_ped_save(:,1),L_ped_save(:,2),L_ped_save(:,3),'.y'),hold on;  %check it
%{
%check points and check the z values
eps = 0.8;
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
        if( L_ped_dx(k,2)< selection_L_body2(j,2)+eps && L_ped_dx(k,2)>= selection_L_body2(j,2)-eps*2 && ...
            L_ped_dx(k,1)< selection_L_body2(j,1)+eps && L_ped_dx(k,1)>= selection_L_body2(j,1)-eps && ...
            L_ped_dx(k,3)< selection_L_body2(j,3)+eps && L_ped_dx(k,3)>= selection_L_body2(j,3)-eps)
         %clean better
        ped_save(count,:) = L_ped_dx(k,:) ; 
        count = count+1;
        
        end 
     end   
end
%plot3(ped_save(:,1),ped_save(:,2),ped_save(:,3),'.y');  %check it
ped_new = setdiff(L_ped_dx,ped_save,'rows');
%plot3(ped_new(:,1), ped_new(:,2), ped_new(:,3),'*r'),hold on;
end
%}
%}
end