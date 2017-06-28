%function to remove the overlap. L_body is the main body and L_ped_dx is
%the body where the points have to be removed

function ped_new =check_overlap_lam( main_body, second_body )       %second body =lam
%remove overlapping points on the pedicles, and on the lamina
%check points and check the z values
% in this case I have to change the condition
eps = 0.5;
count = 1; count_new = 1 ;
if(mean(main_body(:,1)<0) )
    [r,c] = find(second_body(:,1)<=max(main_body(:,1)));
else
    [r,c] = find(second_body(:,1)>=min(main_body(:,1)));
end

selection_L_body= second_body(r,:);
%plot3(selection_L_body(:,1),selection_L_body(:,2),selection_L_body(:,3),'Ob'), hold on, grid on;
clear r; clear c;
[r,c] = find(selection_L_body(:,3)>=min(main_body(:,3)));
selection_L_body2 = selection_L_body(r,:);
%plot3(main_body(:,1),main_body(:,2),main_body(:,3),'.r'), hold on;
%plot3(selection_L_body2(:,1),selection_L_body2(:,2),selection_L_body2(:,3),'Ob'), hold on, grid on;

%delete for difference
ped_new = setdiff(second_body,selection_L_body2,'rows');
%plot3(ped_new(:,1), ped_new(:,2), ped_new(:,3),'*r'),hold on;

end