function [L_body, L_proc, L_ped, L_lam2,EPsup,EPinf] = position_posterior_el(L_body, L_lam, L_proc, EPsup,EPinf,...
    EPWu_half, EPDu, hL, lam_l, sc_d, sc_w, PD_w, PD_h, TP_wu, TP_wi,  PDt, PDs)%F_h,
%{
function [L_body, L_proc, L_ped, L_lam2,vh_vert_regist] = position_posterior_el(L_body, L_lam, L_proc, L_ped,...
    EPWu_half, EPDu, hL, lam_l, sc_d, sc_w, PD_w, PD_h, TP_wu, TP_wi, vh_vert);%F_h,
%}
% function to put in the right position the pedicles, lamina and body.
% Evaluate the min y for the body and the max y for the lamina, this is the
% spinal canal depth. then evaluate the distance between the pedicles, this
% is the spinal canal width
% 
% Do also the scaling for the lamina and lateral processes.

%% -------- v BODY
%selection of points x=0 for the body;
evaluate_body = L_body(L_body(:,1)>=-0.05 & L_body(:,1)<=0.05,: );
       
%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.b'),hold on;
y_min_body = min(evaluate_body(:,2)); % y_min_body is one of the points to evaluate the spinal canal depth
y_max_body = max(evaluate_body(:,2));
%plot3(evaluate_body(:,1),evaluate_body(:,2),evaluate_body(:,3),'*r'),hold on;
%scaling EPDu
init_EPDu = (y_max_body - y_min_body);
sc_EPDu = EPDu/init_EPDu;
L_body = [L_body(:,1),L_body(:,2)*sc_EPDu,L_body(:,3)];
EPsup = [EPsup(:,1),EPsup(:,2)*sc_EPDu,EPsup(:,3)];
EPinf = [EPinf(:,1),EPinf(:,2)*sc_EPDu,EPinf(:,3)];
%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.r'),hold on;


%scaling EPWu
%consider the points on the upper surface, and evaluate the EPW there
evaluate_body = L_body(find(L_body(:,3)>0),:);

%plot3(evaluate_body(:,1),evaluate_body(:,2),evaluate_body(:,3),'.b'),hold on;
x_min_body = min(evaluate_body(:,1)); % x_min_body is one of the points to evaluate the spinal canal depth
x_max_body = max(evaluate_body(:,1));
%plot3(evaluate_body(:,1),evaluate_body(:,2),evaluate_body(:,3),'*r'),hold on;

init_EPWu = (x_max_body - x_min_body);
sc_EPWu = (EPWu_half*2)/init_EPWu;   % evaluate directly EPWu
L_body = [L_body(:,1)*sc_EPWu, L_body(:,2),L_body(:,3)];
EPsup = [EPsup(:,1)*sc_EPWu,EPsup(:,2),EPsup(:,3)];
EPinf = [EPinf(:,1)*sc_EPWu,EPinf(:,2),EPinf(:,3)];
%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.r'),grid on,hold on;

%% ------- LAMINA

evaluate_lam = L_lam(L_lam(:,1)>=-0.5 & L_lam(:,1)<=0.5,:);
       
%plot3(evaluate_lam(:,1),evaluate_lam(:,2),evaluate_lam(:,3),'*b');
y_min_lam = min(evaluate_lam(:,2)); % y_min_lam 
y_max_lam = max(evaluate_lam(:,2)); % y_max_lam is the other point to evaluate the spinal canal depth
init_lam = y_max_lam - y_min_lam ;

sc_lam_y = lam_l/init_lam ;%scaling factor
L_lam =  [ L_lam(:,1), L_lam(:,2)*sc_lam_y, L_lam(:,3) ];


%evaluate again the position and then move the lamina according to the
%canal depth parameter sc_d
evaluate_lam = L_lam(L_lam(:,1)>=-0.5 & L_lam(:,1)<=0.5,:);

%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.g'), hold on;
%plot3(L_body(:,1),L_body(:,2),L_body(:,3),'.g'), hold on;


y_min_lam = min(evaluate_lam(:,2)); % y_min_lam 
y_max_lam = max(evaluate_lam(:,2)); % y_max_lam is the other point to evaluate the spinal canal depth
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.r'),hold on;
mean_body = L_body(L_body(:,1)<0.05 & L_body(:,1)>-0.05,:);
y_min_body = min(mean_body(:,2));
%plot3(0,y_min_body,0,'Ob'),hold on;
sc_d = sc_d;
init_sc_d = y_min_body - y_max_lam;
t_y = sc_d - abs(init_sc_d) ; %define the t value to traslate the lamina
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.r'), hold on;
L_lam = [ L_lam(:,1), L_lam(:,2)+t_y, L_lam(:,3) ]; % the lamina is defined in the x<0 range.
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.r'), hold on;
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.b'),hold on;

clear t_y;

%% --- inferior processes
% inferior processes (belonging to the lamina in this model)
%to select only the inferior part, cut the part for z<0 and x>7 and or x<-7

%find the coordinates of the inferior points of the lamina those belonging
%to the facets
%lam_inf_fac_min = L_lam( L_lam(:,2)<-45 & L_lam(:,2)>-55 & L_lam(:,2)<0,: );
lam_inf_fac_min =  L_lam(L_lam(:,2)<max(L_lam(:,2))-5 & L_lam(:,2)>max(L_lam(:,2))-20 ... 
    & L_lam(:,3)<-5 & L_lam(:,1)<0,: ) ;
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.r'),hold on;
%plot3(lam_inf_fac_min(:,1),lam_inf_fac_min(:,2),lam_inf_fac_min(:,3),'Ob'),hold on;
[minim_min,ind] = min(lam_inf_fac_min(:,1));
%plot3(lam_inf_fac_min(ind,1),lam_inf_fac_min(ind,2),lam_inf_fac_min(ind,3),'Og'),hold on;
dist_TP_wi = abs(minim_min(1)); %symmetric to zero, x distance
sc_TP_wi = (TP_wi/2)/dist_TP_wi; %scaling factor


% do the same procedure for inferior_proc_min and max
cut_min = -5;
cut_max = 5;
z_thres = 0.2;
inferior_proc_max = L_lam(L_lam(:,3)<z_thres & L_lam(:,1)>cut_max,: );
inferior_proc_min = L_lam(L_lam(:,3)<z_thres & L_lam(:,1)<cut_min,: );

L_lam = setdiff(L_lam, inferior_proc_min,'rows'); %remove the part in common
L_lam = setdiff(L_lam, inferior_proc_max,'rows'); %remove the part in common


L_lam = [L_lam(:,1)*1.2, L_lam(:,2), L_lam(:,3)];
%plot3(L_lam(:,1),L_lam(:,2),L_lam(:,3),'.g'), hold on;

%find the min and max x in corrispondence of the central part (y) of the lamina
eps = 0.5;
L_lam_mean = L_lam(L_lam(:,3)<=mean(L_lam(:,3))+eps & L_lam(:,3)>=mean(L_lam(:,3))-eps,:);
cut_min = min(L_lam_mean(:,1));
cut_max = max(L_lam_mean(:,1));

%find a scaling factor for the facets' length
L_body_select = L_body(L_body(:,1)<0.5 & L_body(:,1)>-0.5,: );
z_max = max(L_body_select(:,3));

z_min = min(inferior_proc_max(:,3));
init_dist_f = z_max - z_min;

%scal_z = F_h/init_dist_f;    %scaling factor for all the length; to apply this scaling I have to slightly change the geometry
                            % to do in the next version
%use with the next version
%inferior_proc_min = [inferior_proc_min(:,1)*sc_TP_wi,inferior_proc_min(:,2),inferior_proc_min(:,3)*scal_z];%to do:evaluate this value from abaqus
%inferior_proc_max = [inferior_proc_max(:,1)*sc_TP_wi,inferior_proc_max(:,2),inferior_proc_max(:,3)*scal_z];

inferior_proc_min = [inferior_proc_min(:,1)*sc_TP_wi,inferior_proc_min(:,2),inferior_proc_min(:,3)*1.2];
inferior_proc_max = [inferior_proc_max(:,1)*sc_TP_wi,inferior_proc_max(:,2),inferior_proc_max(:,3)*1.2];


% move back where it has to be
min_pos = max(inferior_proc_min(:,1));
max_pos_z = max(inferior_proc_min(:,3));
%move back in its place
t_x2 = (min_pos) - (cut_min);
%t_z2 = max_pos_z - z_thres;
inferior_proc_min = [inferior_proc_min(:,1)-t_x2, inferior_proc_min(:,2), inferior_proc_min(:,3)];


max_pos = min(inferior_proc_max(:,1));

%move back in its place
t_x2 = (max_pos) - (cut_max);
inferior_proc_max = [inferior_proc_max(:,1)-t_x2, inferior_proc_max(:,2), inferior_proc_min(:,3)];

% remove points in overlap
eps = 0.4;
count =1;
for(j = 1: size(inferior_proc_min))
   for(k = 1:size(L_lam))
    if(L_lam(k,3)<z_thres & inferior_proc_min(j,1)<L_lam(k,1)+eps & ...
            inferior_proc_min(j,1)>L_lam(k,1)-eps)
        remove(count,:) = inferior_proc_min(j,:);
        count = count+1;
    end
   end
end
inferior_proc_min = setdiff(inferior_proc_min,remove,'rows');

eps = 0.4;
count =1;
for(j = 1: size(inferior_proc_max))
   for(k = 1:size(L_lam))
    if(L_lam(k,3)<z_thres & inferior_proc_max(j,1)<L_lam(k,1)+eps & ...
            inferior_proc_max(j,1)>L_lam(k,1)-eps)
        remove(count,:) = inferior_proc_max(j,:);
        count = count+1;
    end
   end
end
inferior_proc_max = setdiff(inferior_proc_max,remove,'rows');


L_lam2 = [L_lam; inferior_proc_min; inferior_proc_max];
%figure
%plot3(L_lam2(:,1), L_lam2(:,2), L_lam2(:,3),'.g'),grid on,hold on;

%% -------------- PEDICLES

% I have to check the positions of the pedicles and adjust them to the
% sc_w. to do it, evaluate the distance of the mean point of each pedicle,
% at the level: y = y_min_body + sc_d/2

mean_body = L_body(L_body(:,1)<0.05 & L_body(:,1)>-0.05,:);
y_min_body = min(mean_body(:,2));
z_max_body = max(mean_body(:,3));


mean_sc = y_min_body - sc_d/2;

P_sag = [ -sc_w/2-PD_w/2 mean_sc]; %x,y
P_trav = [-sc_w/2-PD_w/2 mean_sc z_max_body-1-PD_h/2]; %x y z

angl_sagit = -PDs; 
angl_trasv = 0;%PDt;

L_ped = ellipse_gui(P_sag, P_trav , sc_d,sc_w, PD_h, PD_w, angl_sagit, angl_trasv);


%----------------------------------------

%}
%% ----- PROCESSES

proc_max = L_proc(find(L_proc(:,1)>=0),:);
proc_min = L_proc(find(L_proc(:,1)<=0),:);
%plot3(L_proc(:,1),L_proc(:,2),L_proc(:,3),'.r'),hold on;

%to do: traslate at the same level of the lamina, scale according to TP_wu
max_x = max(proc_max(:,1)); % it s simmetric for x =0; it is the initial TP_wu
min_x = min(proc_max(:,1));

mean_posi_max = mean(proc_max);
mean_posi_min = mean(proc_min);

sc_TP_wu = (TP_wu/2)/max_x;
proc_max = [proc_max(:,1)*sc_TP_wu,proc_max(:,2),proc_max(:,3)];
proc_min =  [proc_min(:,1)*sc_TP_wu,proc_min(:,2),proc_min(:,3)];

mean_posi_max2 = mean(proc_max);
mean_posi_min2 = mean(proc_min);

t_min = mean_posi_min2 - mean_posi_min;
t_max = mean_posi_max2 - mean_posi_max;

proc_max = [proc_max(:,1)+abs(t_max(1)) proc_max(:,2:3)];
proc_min = [proc_min(:,1)-t_min(1) proc_min(:,2:3)];
%plot_matrix(proc_max,'.r');
%plot_matrix(proc_min,'.g');

% find the value to traslate: ealuate the y of the lamina and translate
% there
eps = 0.5;
tmp_lam = L_lam2(find(L_lam2(:,1)>min_x-eps & L_lam2(:,1)<min_x+eps),:);
y_tmp_lam = max(tmp_lam(:,2));
%plot3(L_lam2(:,1),L_lam2(:,2),L_lam2(:,3),'.r'),hold on;
tmp_proc = L_proc(find(L_proc(:,1)>min_x-eps & L_proc(:,1)<min_x+eps),:);
y_tmp_proc = max( tmp_proc(:,2) );

t_y_p = abs(y_tmp_proc - y_tmp_lam);

proc_max = [ proc_max(:,1), proc_max(:,2)-t_y_p, proc_max(:,3)];
proc_min = [ proc_min(:,1), proc_min(:,2)-t_y_p, proc_min(:,3)];
%L_proc = [ proc(:,1), proc(:,2)-t_y_p, proc(:,3)];
L_proc = [proc_min; proc_max];


%check the orientation of the processes
%[L_proc,vh_vert_regist] = orientation_processes_gui(L_proc, L_body,
%vh_vert); % change it to inport the right angle



[L_body, L_ped, L_proc, L_lam2] = double_check_positions_el_gui(L_body, L_ped, L_proc, L_lam2, sc_d);






