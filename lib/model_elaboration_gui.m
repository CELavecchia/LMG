%{
        SCRIPT TO CALL ALL THE FUNCTIONS TOGETHER
- parameterization.m
- load_IVD_2.m
%}
function [mesh_struct_IVD2,L_no] =model_elaboration_gui(dimensions)

%fprintf('fitting\n');

VP_anulus = 0.50;
VP = 50;

%[ L, mesh_struct_IVD ] = parameterization_gui(VP_anulus,dimensions);
fprintf('access the models\n\n')

% Set folder and file name
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','input');
vertName=fullfile(pathName,'L.mat');
ivdName=fullfile(pathName,'mesh_struct_IVD.mat');

outputpathName=fullfile(defaultFolder,'data','output','mat');
outputName=fullfile(outputpathName,'mesh_struct.mat');
cloudpathName=fullfile(defaultFolder,'data','output','point_cloud');

load(vertName);
load(ivdName);

L_body_nt_r = outliers_rem_gui(L.vert_nt.L_body_nt);

% dimensions 
EPWu_half=dimensions.EPWu_half ;
EPDu=dimensions.EPDu; 
EPWi_half=dimensions.EPWi;
EPDi= dimensions.EPDi;
hL=dimensions.hL; 
PDH=dimensions.PDH; 
PDW=dimensions.PDW;
IVD=dimensions.IVD;
lam_l = dimensions.lam_l; 
sc_d=dimensions.sc_d;
sc_w=dimensions.sc_w;
TP_wu=dimensions.TP_wu;
TP_wi=dimensions.TP_wi;
PDt=dimensions.PDt;
PDs=dimensions.PDs;


%% check overlapping points
%plot3(L.vert_nt.L_lam_nt(:,1),L.vert_nt.L_lam_nt(:,2),L.vert_nt.L_lam_nt(:,3),'.r'),hold on;
fprintf('parameterization vertebrae\n\n')
        
ind = 1;
for (j = 1:5)

   [ L_body_rem, L_proc_rem, L_ped_rem,L_lam_rem, EPsup, EPinf] = ...
    remove_overlap_gui(L_body_nt_r(:,ind:ind+2), L.vert_nt.L_proc_grid_nt(:,ind:ind+2),...
    L.vert_nt.L_lam_nt(:,ind:ind+2),L.Vsurf_sup.EP_sup(:,ind:ind+2),L.Vsurf_inf.EP_inf(:,ind:ind+2), ...
    EPWu_half(j), EPDu(j), hL(j), lam_l(j), sc_d(j), sc_w(j), PDW(j), PDH(j),...
    TP_wu(j), TP_wi(j), PDt(j),PDs(j));

    L_no(j).body = L_body_rem;
    L_no(j).proc = L_proc_rem;
    L_no(j).ped = L_ped_rem;
    L_no(j).lam = L_lam_rem;
    L_no(j).EPsup = EPsup;
    L_no(j).EPinf = EPinf;

    ind = ind +3;
    
    %plot3(L_body_rem(:,1),L_body_rem(:,2),L_body_rem(:,3),'.r'),hold on;
    %   plot3(L_proc_rem(:,1),L_proc_rem(:,2),L_proc_rem(:,3),'.r'),hold on;
    %plot3(L_ped_rem(:,1),L_ped_rem(:,2),L_ped_rem(:,3),'.r'),hold on;
    
end

     



%save points/surf 

fprintf('parameterization IVD\n\n')

for(j=1:4)
  %  figure(j),grid on;
   mesh_struct_IVD2(j) = mesh_struct_IVD;
    %[anulus_nodes, nucleus_nodes] = disc_paramet( Disc.anulus(1).V,  Disc.nucleus(1).V, EPWu_half(j), EPDu(j), IVD(j) );
    [mesh_struct_IVD2(j).V] = disc_paramet_new_gui( mesh_struct_IVD.V, EPWu_half(j), EPDu(j), IVD(j) );
    
    %plot_matrix(anulus_nodes,'.');
   
    %save_febio_ivd(Disc.a(j), Disc.n(j), Disc.a(j).fib, Disc.a(j).Eind, j);
   %name = sprintf('./OUTPUT/VP/IVD_VP4_CORR_%i',j);
   
    %save(name,'./output/mesh_struct_IVD2');

    %}  

    %% WRITE STL
    %{
[Ft,Vt]=quad2tri(mesh_struct_IVD2(j).Fb,mesh_struct_IVD2(j).V);
name = sprintf('IVD%d',j)
stlStruct.solidNames={name};
stlStruct.solidVertices={Vt};
stlStruct.solidFaces={Ft};
stlStruct.solidNormals={[]};
fileName=fullfile(sprintf('./ivd%d_MAY_CORR.stl',j));

export_STL_txt(fileName,stlStruct);
%}
    
end
fprintf('saving IVD structure\n\n')

save(outputName,'mesh_struct_IVD2');

%------------------     SAVE  

fprintf('saving the vertebrae point clouds\n');
save_coord_gui( L_no, cloudpathName);
