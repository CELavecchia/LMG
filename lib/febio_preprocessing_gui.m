%% FEBio preprocessing - FU L1-L2-L3
function FEB_struct = febio_preprocessing_gui(modelName, nbodies, Lmes, mesh_struct_IVD2)

%%
% fibres embedded

%clear; close all; clc;
%% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath));

tolerance_par = [0.1];% 0.1 0.5 0.7 0.9];
%contactPenalty = [5];
%vp = [44 50 60];% 44 46 48 50 52 54 56 58 60];
%for(j = 1:1)%length(vp))

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
outputpathName=fullfile(defaultFolder,'data','output','mat');
struct_name = fullfile(outputpathName,'model.mat');

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=20; 

figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=20;
%material properties


%% BONE
Ecort = 12000;
vcort = 0.3;

Ecanc = 100;
vcanc=0.2;

%% CARTILAG eNDPLATE
Ecep = 23.8;
vcep=0.4;

%% NUCLEUS

% LINEAR ELASTIC
%Enucl=1;
%vnucl=0.4;

% Mooney-Rivlin
%c10_nu =0.12; c01_nu =0.03; %Schmidt "Application of a calibration method..."

% Isotropic elastic (Momeni et al)
Enucl = 2; %MPa
vnucl = 0.499; 



%% ANULUS
% isotropic elastic
Eanul=5;
vanul=0.3;


%c10_an =0.18; c01_an =0.045;%Schmidt "Application of a calibration method..."


% Mooney-Rivlin Fan et al
%{
k = 1; c1= 0.56; c2= 0.14;
ksi =1;
%}

% MH Holmes Mow
%{
E_an =0.0233;
v=0.15
beta=2.05;
%}

% hyperelastic HGO (Momeni et al)

c10_an = 0.85; %[MPa]
c01_an = 0;
k = 0; %incompressible term
%fibers coeff Biazial values (O connell, from Momeri et al
ksi = 2.8; %[MPa]
alpha = 90;
beta = 2;
theta =30;
phi=90;
%'ksi','alpha','beta','theta','phi'

%Material parameters


% contact properties
contactPenalty=11;
tolerance = tolerance_par(1);% 0.5;
%{
contactPenalty_up=7;
tolerance_up = 0.4;

contactPenalty_inf=5;
tolerance_inf = 0.6;
%}
%% import datasets
%{
load('C:/Users/lavecchc/Dropbox/Internship_Melbourne/model_November/test_run/FEB_struct_IVD.mat');
load('./struct_ID.mat');         %struct_ID.nodes; struct_ID.ID (all the nodes)
%load('./OUTPUT/struct_layers2.mat');     %in each struct layers layers_struct(j) there 
                               %are the ID arranged in vertical way,from the top to the bottom 
load('./struct_layers2.mat'); 

[supe, infe] = select_nodes_BC(struct_IVD,layers_struct);
%}
%% 
%==============================================================
fprintf('import vertebrae dataset\n\n');

%[Lmes] = FU_whole_model;
%load('C:/Users/lavecchc/Dropbox/Internship_Melbourne/model_November/febio_preprocessing/Lmes_v2.mat');
%for(j=1:nbodies)
%    name = sprintf('./Lmes_v2.mat');
% load(name);
 %Lmes2(j) = Lmes;
%end
%name =sprintf('C:/Users/lavecchc/Dropbox/Internship_Melbourne/model_November/SENS_ANALYSIS/L1L5_VP%d.mat',vp(j)); %IVD
%load(name);


%[mesh_struct_IVD,Lmes] = vertical_orientation(mesh_struct_IVD2,Lmes);

% IVD 1-2
VT = mesh_struct_IVD2(1).V;
E_Ivd1 = mesh_struct_IVD2(1).E;
bC = mesh_struct_IVD2(1).faceBoundaryMarker;
Fbnew = mesh_struct_IVD2(1).Fb;
C = mesh_struct_IVD2(1).C;

% IVD 2-3
VT2 = mesh_struct_IVD2(2).V;
E_Ivd2 = mesh_struct_IVD2(2).E;
bC2 = mesh_struct_IVD2(2).faceBoundaryMarker;
Fbnew2 = mesh_struct_IVD2(2).Fb;
C2 = mesh_struct_IVD2(2).C;

% IVD 3-4
VT3 = mesh_struct_IVD2(3).V;
E_Ivd3 = mesh_struct_IVD2(3).E;
bC3 = mesh_struct_IVD2(3).faceBoundaryMarker;
Fbnew3 = mesh_struct_IVD2(3).Fb;
C3 = mesh_struct_IVD2(3).C;

% IVD 4-5
VT4 = mesh_struct_IVD2(4).V;
E_Ivd4 = mesh_struct_IVD2(4).E;
bC4 = mesh_struct_IVD2(4).faceBoundaryMarker;
Fbnew4 = mesh_struct_IVD2(4).Fb;
C4 = mesh_struct_IVD2(4).C;


%% define nodeset for the attachment points of the ligaments  ----- move this function!!
%ligam = attachment_points(Lmes);


%% join datasets

%[x_node_IVD1, y_node_IVD1, x_node_IVD2, y_node_IVD2] = nodeset(VT,VT2);

if nbodies ==2              % FU L1 - IVD1 - L2
    V = [Lmes(1).VT; Lmes(2).VT; VT];
    E1 = Lmes(1).E;
    E2 = Lmes(2).E + size(Lmes(1).VT,1);
    Lmes(2).Fb = Lmes(2).Fb + size(Lmes(1).VT,1);
    Fbnew = Fbnew + size(Lmes(1).VT,1) + size(Lmes(2).VT,1) ;%IVD1
    febMatL1= Lmes(1).elementMaterialIndices;
    febMatL1(Lmes(1).elementMaterialIndices==-2)=1;
    febMatL1(Lmes(1).elementMaterialIndices==-3)=2;

    febMatL2= Lmes(2).elementMaterialIndices;
    febMatL2(Lmes(2).elementMaterialIndices==-2)=1;
    febMatL2(Lmes(2).elementMaterialIndices==-3)=2;
    
    E_ivd1 = E_Ivd1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1);%+ size(Lmes(4).VT,1)+ size(Lmes(5).VT,1); 
    E_ivd1 = [E_ivd1(:,5:8) E_ivd1(:,1:4)]; %it was inverted for febio
    E_ivd2 = E_Ivd2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(VT,1);% size(Lmes(4).V,1)+ size(Lmes(5).V,1) + size(VT,1); 
    E_ivd2 = [E_ivd2(:,5:8) E_ivd2(:,1:4)];%it was inverted for febio
    
    
   

elseif nbodies ==3          %  L1 - IVD1 - L2 - IVD2 - L3
    V = [Lmes(1).VT; Lmes(2).VT; Lmes(3).VT; VT; VT2];
    E1 = Lmes(1).E;
    E2 = Lmes(2).E + size(Lmes(1).VT,1);
    E3 = Lmes(3).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1);
    Lmes(2).Fb = Lmes(2).Fb + size(Lmes(1).VT,1);
    Lmes(3).Fb = Lmes(3).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1);
    Fbnew = Fbnew + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) ;%IVD1
    Fbnew2 = Fbnew2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1);%IVD2
    febMatL1= Lmes(1).elementMaterialIndices;
    febMatL1(Lmes(1).elementMaterialIndices==-2)=1;
    febMatL1(Lmes(1).elementMaterialIndices==-3)=2;

    febMatL2= Lmes(2).elementMaterialIndices;
    febMatL2(Lmes(2).elementMaterialIndices==-2)=1;
    febMatL2(Lmes(2).elementMaterialIndices==-3)=2;

    febMatL3= Lmes(3).elementMaterialIndices;
    febMatL3(Lmes(3).elementMaterialIndices==-2)=1;
    febMatL3(Lmes(3).elementMaterialIndices==-3)=2;
    
    E_ivd1 = E_Ivd1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);%+ size(Lmes(4).VT,1)+ size(Lmes(5).VT,1); 
    E_ivd1 = [E_ivd1(:,5:8) E_ivd1(:,1:4)]; %it was inverted for febio
    E_ivd2 = E_Ivd2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(VT,1);% size(Lmes(4).V,1)+ size(Lmes(5).V,1) + size(VT,1); 
    E_ivd2 = [E_ivd2(:,5:8) E_ivd2(:,1:4)];%it was inverted for febio

  
elseif nbodies ==4          %  L1 - IVD1 - L2 - IVD2 - L3 - IVD3 - L4
    V = [Lmes(1).VT; Lmes(2).VT; Lmes(3).VT; Lmes(4).VT; VT; VT2; VT3];
    E1 = Lmes(1).E;
    E2 = Lmes(2).E + size(Lmes(1).VT,1);
    E3 = Lmes(3).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1);
    E4 = Lmes(4).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1) + size(Lmes(3).VT,1);
    Lmes(2).Fb = Lmes(2).Fb + size(Lmes(1).VT,1);
    Lmes(3).Fb = Lmes(3).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1);
    Lmes(4).Fb = Lmes(4).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1)+size(Lmes(3).VT,1);
    Fbnew = Fbnew + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) ;%IVD1
    Fbnew2 = Fbnew2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1);%IVD2
    Fbnew3 = Fbnew3 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1)+ size(VT2,1);%IVD3
    
    febMatL1= Lmes(1).elementMaterialIndices;
    febMatL1(Lmes(1).elementMaterialIndices==-2)=1;
    febMatL1(Lmes(1).elementMaterialIndices==-3)=2;

    febMatL2= Lmes(2).elementMaterialIndices;
    febMatL2(Lmes(2).elementMaterialIndices==-2)=1;
    febMatL2(Lmes(2).elementMaterialIndices==-3)=2;

    febMatL3= Lmes(3).elementMaterialIndices;
    febMatL3(Lmes(3).elementMaterialIndices==-2)=1;
    febMatL3(Lmes(3).elementMaterialIndices==-3)=2;
    
    febMatL4= Lmes(4).elementMaterialIndices;
    febMatL4(Lmes(4).elementMaterialIndices==-2)=1;
    febMatL4(Lmes(4).elementMaterialIndices==-3)=2;
    
    E_ivd1 = E_Ivd1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1); 
    E_ivd1 = [E_ivd1(:,5:8) E_ivd1(:,1:4)]; %it was inverted for febio
    E_ivd2 = E_Ivd2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1)+ size(VT,1)+ size(VT2,1);
    E_ivd2 = [E_ivd2(:,5:8) E_ivd2(:,1:4)];%it was inverted for febio
    E_ivd3 = E_Ivd3 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1)+ size(VT,1)+ size(VT2,1)+ size(VT3,1);
    E_ivd3 = [E_ivd3(:,5:8) E_ivd3(:,1:4)];%it was inverted for febio
    

    
end
    
%nodeset for the output
%{
x_node_IVD1 =  x_node_IVD1 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1);
y_node_IVD1 = y_node_IVD1 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1);
x_node_IVD2 = x_node_IVD2 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1) +size(VT(:,1),1);
y_node_IVD2 = y_node_IVD2 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1) +size(VT(:,1),1);
%}


%{
E1 = Lmes(1).E;
E2 = Lmes(2).E + size(Lmes(1).VT,1);
E3 = Lmes(3).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1);
%}
%E4 = Lmes(4).E1 + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);
%E5 = Lmes(5).E1 + size(Lmes(1).V,1)+ size(Lme2(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1);
%{
E_ivd1 = E_Ivd1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);%+ size(Lmes(4).VT,1)+ size(Lmes(5).VT,1); 
E_ivd1 = [E_ivd1(:,5:8) E_ivd1(:,1:4)]; %it was inverted for febio
E_ivd2 = E_Ivd2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(VT,1);% size(Lmes(4).V,1)+ size(Lmes(5).V,1) + size(VT,1); 
E_ivd2 = [E_ivd2(:,5:8) E_ivd2(:,1:4)];%it was inverted for febio
%}
%E_ivd3 = E_Ivd3 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);%+ size(Lmes(4).V,1)+ size(Lmes(5).V,1) + size(VT,1)+ size(VT2,1); 

%E_ivd3 = [E_ivd3(:,5:8) E_ivd3(:,1:4)];%it was inverted for febio

%E_ivd4 = E_Ivd4 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);%+ size(Lmes(4).VT,1)+ size(Lmes(5).V,1) + size(VT,1)+ size(VT2,1)+ size(VT3,1); 
%E_ivd4 = [E_ivd4(:,5:8) E_ivd4(:,1:4)];%it was inverted for febio

%join faces as well
%{
Lmes(2).Fb = Lmes(2).Fb + size(Lmes(1).VT,1);
Lmes(3).Fb = Lmes(3).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1);
%}
%Lmes(4).Fb = Lmes(4).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1)+ size(Lmes(3).VT,1);
%Lmes(5).Fb = Lmes(5).Fb + size(Lmes(1).V,1)+size(Lmes(2).V,1)+ size(Lmes(3).V,1)+size(Lmes(4).V,1);

%{
Fbnew = Fbnew + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) ;%IVD1
Fbnew2 = Fbnew2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1);%IVD2
%Fbnew3 = Fbnew3 + size(Lmes(1).V,1) + size(Lmes(2).V,1)+ size(Lmes(3).V,1) + size(VT2,1);%IVD2
%Fbnew4 = Fbnew4 + size(Lmes(1).V,1) + size(Lmes(2).V,1)+ size(Lmes(3).V,1) + size(VT3,1);%IVD2
%}

% update indices of the ligaments attaching points
%ligam2 = update_indeces(ligam, size(Lmes(1).VT(:,1),1), size(Lmes(2).VT(:,1)) );

%materials indeces

%{
febMatL1= Lmes(1).elementMaterialIndices;
febMatL1(Lmes(1).elementMaterialIndices==-2)=1;
febMatL1(Lmes(1).elementMaterialIndices==-3)=2;

febMatL2= Lmes(2).elementMaterialIndices;
febMatL2(Lmes(2).elementMaterialIndices==-2)=1;
febMatL2(Lmes(2).elementMaterialIndices==-3)=2;

febMatL3= Lmes(3).elementMaterialIndices;
febMatL3(Lmes(3).elementMaterialIndices==-2)=1;
febMatL3(Lmes(3).elementMaterialIndices==-3)=2;
%}
%{
febMatL4= Lmes(4).elementMaterialIndices;
febMatL4(Lmes(4).elementMaterialIndices==-2)=1;
febMatL4(Lmes(4).elementMaterialIndices==-3)=2;

febMatL5= Lmes(5).elementMaterialIndices;
febMatL5(Lmes(5).elementMaterialIndices==-2)=1;
febMatL5(Lmes(5).elementMaterialIndices==-3)=2;
%}



%{
febMatIVD3= C3;
febMatIVD3(C32==1)=3;
febMatIVD3(C3==2)=4;

febMatIVD4= C4;
febMatIVD4(C4==1)=3;
febMatIVD4(C4==2)=4;
%}
%element indices for the IVD
febMatIVD= C;
%anulus   % you can select also different material prop for each layer
febMatIVD(C==2)=3; %nucleus
febMatIVD(C==3)=4;
febMatIVD(C==4)=3;
febMatIVD(C==5)=4;
febMatIVD(C==6)=3;
febMatIVD(C==7)=4;
febMatIVD(C==8)=3;
febMatIVD(C==9)=4;
%nucleus
febMatIVD(C==10)=5;
febMatIVD(C==11)=5;


%IVD2
febMatIVD2= C2;
febMatIVD2(C2==2)=3; %nucleus
febMatIVD2(C2==3)=4;
febMatIVD2(C2==4)=3;
febMatIVD2(C2==5)=4;
febMatIVD2(C2==6)=3;
febMatIVD2(C2==7)=4;
febMatIVD2(C2==8)=3;
febMatIVD2(C2==9)=4;
%nucleus
febMatIVD2(C2==10)=5;
febMatIVD2(C2==11)=5;

% IVD3
febMatIVD3= C3;
febMatIVD3(C3==2)=3; %nucleus
febMatIVD3(C3==3)=4;
febMatIVD3(C3==4)=3;
febMatIVD3(C3==5)=4;
febMatIVD3(C3==6)=3;
febMatIVD3(C3==7)=4;
febMatIVD3(C3==8)=3;
febMatIVD3(C3==9)=4;
%nucleus
febMatIVD3(C3==10)=5;
febMatIVD3(C3==11)=5;


%% CREATING FIBRE DIRECTIONS
%{
%% IVD1
%Compute coordinates at the centre of elements for fibre origings
X=V(:,1); Y=V(:,2); Z=V(:,3);

%elements anulus layer1
E_anulus_l1 = E_ivd1(febMatIVD==3,:);%select only the anulus elements
index_l1 = find(febMatIVD==3);
%index_l1_old = index_l1;
%index_l1 = index_l1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(3).VT,1);
index_l1 = index_l1 + size(Lmes(1).E,1) + size(Lmes(2).E,1)+ size(Lmes(3).E,1) ;


E_nucleus = E_ivd1(febMatIVD==5,:);
febMatNucleus = ones(size(E_nucleus,1),1);

E_save =select_elements_endplates( E_nucleus, bC, Fbnew );

%find the index of these elements in the original E_ivd1 and assign for the
%cartilagineus endplates
[val, index_EIV, index_Esav] = intersect(E_ivd1,E_save,'rows','stable');
febMatIVD(index_EIV)=6;


XE=mean(X(E_anulus_l1),2); YE=mean(Y(E_anulus_l1),2); ZE=mean(Z(E_anulus_l1),2);
VE=[XE(:) YE(:) ZE(:)];

%elementMaterialIndices=febMatIVD;

FI1=3*ones(length(E_anulus_l1),1);%elementMaterialIndices==3; %anulus
%FI2=elementMaterialIndices==5;  %nucleus

zDir=[0 0 1];
VE_XY=VE;
VE_XY(:,3)=0;
VE_XY=vecnormalize(VE_XY);
VF=cross(VE_XY,zDir(ones(size(VE,1),1),:)); %I.e. radial tangent
VF=vecnormalize(VF);
%VF(FI1,:)=-VF(FI1,:); %I.e. radial tangent flipped
%VF(:,3)=1; 
%VF=vecnormalize(VF);

E_anulus_l2 = E_ivd1(febMatIVD==4,:);%select only the anulus elements
index_l2 = find(febMatIVD==4);
%index_l2 = index_l2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(3).VT,1)+size(index_l1,1);
index_l2 = index_l2 + size(Lmes(1).E,1) + size(Lmes(2).E,1)+ size(Lmes(3).E,1);%+size(index_l1,1);


XE=mean(X(E_anulus_l2),2); YE=mean(Y(E_anulus_l2),2); ZE=mean(Z(E_anulus_l2),2);
VE2=[XE(:) YE(:) ZE(:)];

elementMaterialIndices=febMatIVD;

FI12=4*ones(length(E_anulus_l2),1);%elementMaterialIndices==3; %anulus
%FI2=elementMaterialIndices==4;  %nucleus

zDir=[0 0 -1];
VE_XY2=VE2;
VE_XY2(:,3)=0;
VE_XY2=vecnormalize(VE_XY2);
VF2=cross(VE_XY2,zDir(ones(size(VE2,1),1),:)); %I.e. radial tangent
VF2=vecnormalize(VF2);
%VF2(FI12,:)=-VF2(FI12,:); %I.e. radial tangent flipped
%VF2(:,3)=1; 
%VF2=vecnormalize(VF2);


[Ff,Vf,Cf]=quiver3Dpatch(VE(:,1),VE(:,2),VE(:,3),VF(:,1),VF(:,2),VF(:,3),ones(size(VF,1),1),[1 1]);
[Ff2,Vf2,Cf2]=quiver3Dpatch(VE2(:,1),VE2(:,2),VE2(:,3),VF2(:,1),VF2(:,2),VF2(:,3),ones(size(VF2,1),1),[1 1]);
%{
hf1=cFigure;
title('Fibre directions ','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
%hp=patch('Faces',quad,'Vertices',VT(:,1:3),'EdgeColor','none','FaceColor','flat','CData',2,'FaceAlpha',faceAlpha2);
hp=patch('Faces',E_anulus_l2,'Vertices',V(:,1:3),'CData',2,'EdgeColor','b','FaceAlpha',0.1);
hp=patch('Faces',E_anulus_l1,'Vertices',V(:,1:3),'CData',2,'EdgeColor','g','FaceAlpha',0.1);
hp=patch('Faces',E_nucleus,'Vertices',V(:,1:3),'CData',2,'EdgeColor','r','FaceAlpha',0.1);
hp=patch('Faces',E_save,'Vertices',V(:,1:3),'CData',2,'EdgeColor','c','FaceAlpha',0.1);
plotV(VE,'k.');
hp=patch('Faces',Ff,'Vertices',Vf,'EdgeColor','none','FaceColor','r');
hp=patch('Faces',Ff2,'Vertices',Vf2,'EdgeColor','none','FaceColor','r');
%}
%% IVD2
%Compute coordinates at the centre of elements for fibre origings
%X=VT2(:,1); Y=VT2(:,2); Z=VT2(:,3);
%X=V(:,1); Y=VT2(:,2); Z=VT2(:,3);

%elements anulus layer1
E_anulus2_l1 = E_ivd2(febMatIVD2==3,:);%select only the anulus elements
index2_l1 = find(febMatIVD2==3);
%index2_l1 = index2_l1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(3).VT,1)+ size(E_ivd1,1);
index2_l1 = index2_l1 + size(Lmes(1).E,1) + size(Lmes(2).E,1)+ size(Lmes(3).E,1)+ size(E_ivd1,1);

E_nucleus2 = E_ivd2(febMatIVD2==5,:);
XE=mean(X(E_anulus2_l1),2); YE=mean(Y(E_anulus2_l1),2); ZE=mean(Z(E_anulus2_l1),2);
VE3=[XE(:) YE(:) ZE(:)];


E_save2 =select_elements_endplates( E_nucleus2, bC2, Fbnew2 );

%find the index of these elements in the original E_ivd1
[val, index_EIV2, index_Esav2] = intersect(E_ivd2,E_save2,'rows','stable');
febMatIVD2(index_EIV2)=6;

elementMaterialIndices=febMatIVD;

FI3=3*ones(length(E_anulus_l1),1);%elementMaterialIndices==3; %anulus
FI=elementMaterialIndices==5;  %nucleus

%select_elements_endplates( E_nucleus, febMatIVD, VT, bC );

zDir=[0 0 1];
VE_XY3=VE3;
VE_XY3(:,3)=0;
VE_XY3=vecnormalize(VE_XY3);
VF3=cross(VE_XY3,zDir(ones(size(VE3,1),1),:)); %I.e. radial tangent
VF3=vecnormalize(VF3);
%{
VF3(FI1,:)=-VF3(FI1,:); %I.e. radial tangent flipped
VF3(:,3)=1; 
VF3=vecnormalize(VF3);
%}

E_anulus2_l2 = E_ivd2(febMatIVD2==4,:);%select only the anulus elements
index2_l2 = find(febMatIVD2==4);
%index2_l2 = index2_l2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(3).VT,1)+ size(E_ivd1,1)+size(index2_l1,1);
index2_l2 = index2_l2 + size(Lmes(1).E,1) + size(Lmes(2).E,1)+ size(Lmes(3).E,1)+ size(E_ivd1,1);%+size(index2_l1,1);


%E_nucleus = E_Ivd1(febMatIVD==4,:);
XE=mean(X(E_anulus2_l2),2); YE=mean(Y(E_anulus2_l2),2); ZE=mean(Z(E_anulus2_l2),2);
VE4=[XE(:) YE(:) ZE(:)];

elementMaterialIndices=febMatIVD2;

FI4=4*ones(length(E_anulus2_l2),1);%elementMaterialIndices==3; %anulus
%FI2=elementMaterialIndices==4;  %nucleus

%select_elements_endplates( E_nucleus, febMatIVD, VT, bC );

zDir=[0 0 -1];
VE_XY4=VE4;
VE_XY4(:,3)=0;
VE_XY4=vecnormalize(VE_XY4);
VF4=cross(VE_XY4,zDir(ones(size(VE4,1),1),:)); %I.e. radial tangent
VF4=vecnormalize(VF4);
%{
VF4(FI4,:)=-VF4(FI4,:); %I.e. radial tangent flipped
VF4(:,3)=1; 
VF4=vecnormalize(VF4);
%}
%Create patch data for fibres
[Ff,Vf,Cf]=quiver3Dpatch(VE(:,1),VE(:,2),VE(:,3),VF(:,1),VF(:,2),VF(:,3),ones(size(VF,1),1),[1 1]);
[Ff2,Vf2,Cf2]=quiver3Dpatch(VE2(:,1),VE2(:,2),VE2(:,3),VF2(:,1),VF2(:,2),VF2(:,3),ones(size(VF2,1),1),[1 1]);
[Ff3,Vf3,Cf3]=quiver3Dpatch(VE3(:,1),VE3(:,2),VE3(:,3),VF3(:,1),VF3(:,2),VF3(:,3),ones(size(VF3,1),1),[1 1]);
[Ff4,Vf4,Cf4]=quiver3Dpatch(VE4(:,1),VE4(:,2),VE4(:,3),VF4(:,1),VF4(:,2),VF4(:,3),ones(size(VF4,1),1),[1 1]);



E_ind_vect_fibres = [Vf;Vf2;Vf3;Vf4];
%E_ind_fibres = [E_anulus_l1;E_anulus_l2;E_anulus2_l1;E_anulus2_l1];
E_ind_fibres = [index_l1;index_l2;index2_l1;index2_l2];


%Create patch data for fibres
[Ff,Vf,Cf]=quiver3Dpatch(VE(:,1),VE(:,2),VE(:,3),VF(:,1),VF(:,2),VF(:,3),ones(size(VF,1),1),[1 1]);
[Ff2,Vf2,Cf2]=quiver3Dpatch(VE2(:,1),VE2(:,2),VE2(:,3),VF2(:,1),VF2(:,2),VF2(:,3),ones(size(VF2,1),1),[1 1]);
%{
hf1=cFigure;
title('Fibre directions ','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
%hp=patch('Faces',quad,'Vertices',VT(:,1:3),'EdgeColor','none','FaceColor','flat','CData',2,'FaceAlpha',faceAlpha2);
hp=patch('Faces',E_Ivd1,'Vertices',VT(:,1:3),'CData',2,'EdgeColor','c','FaceAlpha',0.2);

hp=patch('Faces',Ff,'Vertices',Vf,'EdgeColor','none','FaceColor','r');
hp=patch('Faces',Ff2,'Vertices',Vf2,'EdgeColor','none','FaceColor','b');
%hp=patch('Faces',Ff3,'Vertices',Vf3,'EdgeColor','none','FaceColor','g');
%hp=patch('Faces',Ff4,'Vertices',Vf4,'EdgeColor','none','FaceColor','c');
%hp=patch('Faces',Ff2,'Vertices',Vf2,'EdgeColor','none','FaceColor','r');
% plotV(VE,'k.');
axis equal; view(3); axis tight; set(gca,'FontSize',fontSize);
%colormap hsv; colorbar; 
set(gca,'FontSize',fontSize);
drawnow;
%}
%}

%% save the model


model.V = V; %vert
    model.E1 = E1; model.E2 = E2; 
    model.E_ivd1 = E_ivd1; 
    model.Fb1 = Lmes(1).Fb; model.Fb2 = Lmes(2).Fb;
    model.matind1 =febMatL1; model.matind2 =febMatL2;
    model.Cb = Lmes(1).Cb; model.Cb2 = Lmes(2).Cb ; 

    model.FbI1 = febMatIVD; 
    model.Fbnew = Fbnew; model.Fbnew2 = Fbnew2; %faces
    model.bcI1 = mesh_struct_IVD2(1).faceBoundaryMarker; 
    %model.E_ind_vect_fibres = E_ind_vect_fibres; model.E_ind_fibres = E_ind_fibres; 
    
    if nbodies ==3
         
         model.E3 = E3;model.E_ivd2 = E_ivd2; %elements
         model.matind3 =febMatL3; %mat ind vert
         model.FbI2 = febMatIVD2;    %mat ind IVD
         model.bcI2 = mesh_struct_IVD2(2).faceBoundaryMarker; %bc IVD
          model.Fb3 = Lmes(3).Fb; %faces vert
          model.Cb3 = Lmes(3).Cb; % cb vert
    elseif nbodies ==4
         model.E4 = E4;model.E_ivd3 = E_ivd3; %elements
         model.matind4 =febMatL4; %mat ind vert
         model.FbI3 = febMatIVD3;    %mat ind IVD
         model.bcI3 = mesh_struct_IVD2(3).faceBoundaryMarker; %bc IVD
          model.Fb4 = Lmes(4).Fb; %faces vert
          model.Cb4 = Lmes(4).Cb; % cb vert
    elseif nbodies ==5
         model.E5 = E5;model.E_ivd4 = E_ivd4; %elements
         model.matind5 =febMatL5; %mat ind vert
         model.FbI4 = febMatIVD4;    %mat ind IVD
         model.bcI4 = mesh_struct_IVD2(4).faceBoundaryMarker; %bc IVD
          model.Fb5 = Lmes(5).Fb; %faces vert
          model.Cb5 = Lmes(5).Cb; % cb vert 
        
        
    end
    
    save(struct_name,'model');
    
    
      

%%  BC
%IVD
logicFace=bC==2;
FI_sup1=Fbnew(logicFace,:);
disc_sup1=unique(FI_sup1(:));
logicFace=bC==1;
FI_inf1=Fbnew(logicFace,:); 
disc_inf1=unique(FI_inf1(:));

logicFace=bC2==2;
FI_sup2=Fbnew2(logicFace,:);
disc_sup2=unique(FI_sup2(:));
logicFace=bC2==1;
FI_inf2=Fbnew2(logicFace,:); 
disc_inf2=unique(FI_inf2(:));
%{
logicFace=bC3==1;
FI_sup3=Fbnew3(logicFace,:);
disc_sup3=unique(FI_sup3(:));
logicFace=bC3==2;
FI_inf3=Fbnew3(logicFace,:); 
disc_inf3=unique(FI_inf2(:));

logicFace=bC4==1;
FI_sup4=Fbnew4(logicFace,:);
disc_sup4=unique(FI_sup4(:));
logicFace=bC4==2;
FI_inf4=Fbnew4(logicFace,:); 
disc_inf4=unique(FI_inf4(:));
%}

%{
logicFace=faceBoundaryMarkerIVD==6;
FI_sup3=Fbnew(logicFace,:);
disc_sup3=unique(FI_sup1(:));
logicFace=faceBoundaryMarkerIVD==5;
FI_inf3=Fbnew(logicFace,:); 
disc_inf3=unique(FI_inf1(:));

logicFace=faceBoundaryMarkerIVD2==6;
FI_sup4=Fbnew2(logicFace,:);
disc_sup4=unique(FI_sup2(:));
logicFace=faceBoundaryMarkerIVD2==5;
FI_inf2=Fbnew2(logicFace,:); 
disc_inf2=unique(FI_inf2(:));
%}
%vertebrae
sz = 0;
for(j =1:nbodies)
    logicFace=Lmes(j).Cb==2; % select only the face outside
    Face_ext= Lmes(j).Fb(logicFace,:);

    [N]=patchNormal(Face_ext, V); 
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
    logicFace2=faceBoundaryMarker==6;
    F(j).vert_inf=Face_ext(logicFace2,:);
    bc(j).vert_inf=unique(F(j).vert_inf(:));

    
sz = size(Lmes(j).VT,1);
%hold on;
%patch('Faces',F_selec2,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);

%patch('Faces',Lmes(j).Fb,'Vertices',V,'FaceColor','flat','CData',Lmes(j).Cb,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
%hold on;
%patch('Faces',Face_ext,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);

%plotV(V(bc(j).vert_sup,:),'b.','MarkerSize',markerSize);hold on;
%plotV(V(bc(j).vert_inf,:),'r.','MarkerSize',markerSize);hold on;
    
end

%% contacts
%vertebrae
%nodeset


%IVD

%I1
logicContactSurf_sup_IVD1 = disc_sup1;
logicContactSurf_inf_IVD1 = disc_inf1;
%I2
logicContactSurf_sup_IVD2 = disc_sup2;
logicContactSurf_inf_IVD2 = disc_inf2;
%I2
%{
logicContactSurf_sup_IVD3 = disc_sup3;
logicContactSurf_inf_IVD3 = disc_inf3;
%I2
logicContactSurf_sup_IVD4 = disc_sup4;
logicContactSurf_inf_IVD4 = disc_inf4;
%}

%% CONSTRUCTING FEB MODEL
  FEB_struct.febio_spec.version='2.0';
    FEB_struct.Module.Type='solid';

    % Defining file names
    FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
    FEB_struct.run_logname=[modelName,'.log']; %FEBio log file name

    %Geometry section
    FEB_struct.Geometry.Nodes=V;
    if nbodies ==2
        FEB_struct.Geometry.Elements={E1 E2 E_ivd1}; %The element sets
        FEB_struct.Geometry.ElementType={'tet4','tet4','hex8'}; %The element types
        FEB_struct.Geometry.ElementMat={febMatL1', febMatL2', febMatIVD'};
        FEB_struct.Geometry.ElementsPartName={'L1' 'L2' 'IVD1'};

    elseif nbodies==3
            FEB_struct.Geometry.Elements={E1 E2 E3 E_ivd1 E_ivd2}; %The element sets
        FEB_struct.Geometry.ElementType={'tet4','tet4','tet4','hex8','hex8'}; %The element types
        FEB_struct.Geometry.ElementMat={febMatL1', febMatL2',febMatL3', febMatIVD', febMatIVD2'};
        FEB_struct.Geometry.ElementsPartName={'L1' 'L2' 'L3' 'IVD1' 'IVD2'};
    end

%% Material section

FEB_struct.Materials{1}.Type = 'isotropic elastic';
FEB_struct.Materials{1}.Name = 'bone';
FEB_struct.Materials{1}.Properties = {'E','v'};
FEB_struct.Materials{1}.Values = {Ecort,vcort};

FEB_struct.Materials{2}.Type = 'isotropic elastic';
FEB_struct.Materials{2}.Name = 'canc';
FEB_struct.Materials{2}.Properties = {'E','v'};
FEB_struct.Materials{2}.Values = {Ecanc,vcanc};


FEB_struct.Materials{3}.Type='solid mixture';
%{
FEB_struct.Materials{3}.Solid{1}.Type='Holmes-Mow';%'Ogden unconstrained';
FEB_struct.Materials{3}.Solid{1}.Properties={'E','v','beta'};%{'c1','m1','cp'};
FEB_struct.Materials{3}.Solid{1}.Values= {E_an,v,beta};%{c1_1,m1_1,k_1};
%}

FEB_struct.Materials{3}.Solid{1}.Type='neo-Hookean';%'Ogden unconstrained';
FEB_struct.Materials{3}.Solid{1}.Properties={'E','v'};%{'c1','m1','cp'};
FEB_struct.Materials{3}.Solid{1}.Values= {Eanul,vanul};%{c1_1,m1_1,k_1};
%{
FEB_struct.Materials{3}.Solid{1}.Type='Mooney-Rivlin';%'Ogden unconstrained';
FEB_struct.Materials{3}.Solid{1}.Properties={'c1','c2','k'};%{'c1','m1','cp'};
FEB_struct.Materials{3}.Solid{1}.Values= {c10_an,c01_an,0};%{c1_1,m1_1,k_1};
%}
 FEB_struct.Materials{3}.Solid{2}.Type='fiber-exp-pow';
 FEB_struct.Materials{3}.Solid{2}.Properties={'ksi','alpha','beta','theta','phi'};
 FEB_struct.Materials{3}.Solid{2}.Values={ksi,alpha,beta,theta,phi};
 FEB_struct.Materials{3}.Solid{2}.AnisoType='mat_axis';
 
  FEB_struct.Materials{3}.Solid{3}.Type='fiber-exp-pow';
 FEB_struct.Materials{3}.Solid{3}.Properties={'ksi','alpha','beta','theta','phi'};
 FEB_struct.Materials{3}.Solid{3}.Values={ksi,alpha,beta,-theta,phi};
 FEB_struct.Materials{3}.Solid{3}.AnisoType='mat_axis';
  
% FEB_struct.Materials{4} = FEB_struct.Materials{3};
 
FEB_struct.Materials{4}.Type='solid mixture';

%{
FEB_struct.Materials{4}.Solid{1}.Type='Holmes-Mow';%'Ogden unconstrained';
FEB_struct.Materials{4}.Solid{1}.Properties={'E','v','beta'};%{'c1','m1','cp'};
FEB_struct.Materials{4}.Solid{1}.Values= {E_an,v,beta};%{c1_1,m1_1,k_1};
%}
FEB_struct.Materials{4}.Solid{1}.Type='neo-Hookean';%'Ogden unconstrained';
FEB_struct.Materials{4}.Solid{1}.Properties={'E','v'};%{'c1','m1','cp'};
FEB_struct.Materials{4}.Solid{1}.Values= {Eanul,vanul};%{c1_1,m1_1,k_1};
%{
FEB_struct.Materials{4}.Solid{1}.Type='Mooney-Rivlin';%'Ogden unconstrained';
FEB_struct.Materials{4}.Solid{1}.Properties={'c1','c2','k'};%{'c1','m1','cp'};
FEB_struct.Materials{4}.Solid{1}.Values= {c10_an,c01_an,0};%{c1_1,m1_1,k_1};
%}
FEB_struct.Materials{4}.Solid{2}.Type='fiber-exp-pow';
 FEB_struct.Materials{4}.Solid{2}.Properties={'ksi','alpha','beta','theta','phi'};
 FEB_struct.Materials{4}.Solid{2}.Values={ksi,alpha,beta,theta,phi};
 FEB_struct.Materials{4}.Solid{2}.AnisoType='mat_axis';
 
  FEB_struct.Materials{4}.Solid{3}.Type='fiber-exp-pow';
 FEB_struct.Materials{4}.Solid{3}.Properties={'ksi','alpha','beta','theta','phi'};
 FEB_struct.Materials{4}.Solid{3}.Values={ksi,alpha,beta,-theta,phi};
 FEB_struct.Materials{4}.Solid{3}.AnisoType='mat_axis';
 
 
 
 %nucleus
 %isotropic linear
 FEB_struct.Materials{5}.Type = 'isotropic elastic';
FEB_struct.Materials{5}.Name = 'nucleus'; %np1
FEB_struct.Materials{5}.Properties = {'E','v'};
FEB_struct.Materials{5}.Values = {Enucl,vnucl};
 
 %{
FEB_struct.Materials{5}.Type = 'neo-Hookean';
FEB_struct.Materials{5}.Name = 'nucleus'; %np1
FEB_struct.Materials{5}.Properties = {'E','v'};
FEB_struct.Materials{5}.Values = {Enucl,vnucl};
 %}
 
 %HM same properties of the anulus
 %{
FEB_struct.Materials{5}.Type='Holmes-Mow';%'Ogden unconstrained';
FEB_struct.Materials{5}.Properties={'E','v','beta'};%{'c1','m1','cp'};
FEB_struct.Materials{5}.Values= {E_an,v,beta};%{c1_1,m1_1,k_1};
%}
%{
FEB_struct.Materials{5}.Type = 'Mooney-Rivlin';
FEB_struct.Materials{5}.Name = 'nucleus'; %np1
FEB_struct.Materials{5}.Properties = {'c1','c2','k'};
FEB_struct.Materials{5}.Values = {c10_nu,c01_nu,0};
%}
 %CEP
FEB_struct.Materials{6}.Type = 'isotropic elastic';
FEB_struct.Materials{6}.Name = 'bone';
FEB_struct.Materials{6}.Properties = {'E','v'};
FEB_struct.Materials{6}.Values = {Ecep,vcep};
%% FIBRES DEFINITION ---> RE-CHECK!



%% BCs

%Defining node sets
%IVD1
FEB_struct.Geometry.NodeSet{1}.Set = disc_sup1; 
FEB_struct.Geometry.NodeSet{1}.Name = 'disc_sup1';
FEB_struct.Geometry.NodeSet{2}.Set = disc_inf1; 
FEB_struct.Geometry.NodeSet{2}.Name = 'disc_inf1';



if nbodies==3
    %IVD2
    FEB_struct.Geometry.NodeSet{3}.Set = disc_sup2; 
    FEB_struct.Geometry.NodeSet{3}.Name = 'disc_sup2';
    FEB_struct.Geometry.NodeSet{4}.Set = disc_inf2; 
    FEB_struct.Geometry.NodeSet{4}.Name = 'disc_inf2';


    %L1
    FEB_struct.Geometry.NodeSet{5}.Set = bc(1).vert_sup; 
    FEB_struct.Geometry.NodeSet{5}.Name = 'vert_L1_sup';
    FEB_struct.Geometry.NodeSet{6}.Set = bc(1).vert_inf; 
    FEB_struct.Geometry.NodeSet{6}.Name = 'vert_L1_inf';
    %L2
    FEB_struct.Geometry.NodeSet{7}.Set = bc(2).vert_sup; 
    FEB_struct.Geometry.NodeSet{7}.Name = 'vert_L2_sup';
    FEB_struct.Geometry.NodeSet{8}.Set = bc(2).vert_inf; 
    FEB_struct.Geometry.NodeSet{8}.Name = 'vert_L2_inf';
    %L3

         FEB_struct.Geometry.NodeSet{9}.Set = bc(3).vert_sup; 
        FEB_struct.Geometry.NodeSet{9}.Name = 'vert_L3_sup';
        FEB_struct.Geometry.NodeSet{10}.Set = bc(3).vert_inf; 
        FEB_struct.Geometry.NodeSet{10}.Name = 'vert_L3_inf';
end

 %L1
    FEB_struct.Geometry.NodeSet{3}.Set = bc(1).vert_sup; 
    FEB_struct.Geometry.NodeSet{3}.Name = 'vert_L1_sup';
    FEB_struct.Geometry.NodeSet{4}.Set = bc(1).vert_inf; 
    FEB_struct.Geometry.NodeSet{4}.Name = 'vert_L1_inf';
    %L2
    FEB_struct.Geometry.NodeSet{5}.Set = bc(2).vert_sup; 
    FEB_struct.Geometry.NodeSet{5}.Name = 'vert_L2_sup';
    FEB_struct.Geometry.NodeSet{6}.Set = bc(2).vert_inf; 
    FEB_struct.Geometry.NodeSet{6}.Name = 'vert_L2_inf';
    %L3

%{
%L4
FEB_struct.Geometry.NodeSet{15}.Set = bc(4).vert_sup; 
FEB_struct.Geometry.NodeSet{15}.Name = 'vert_L4_sup';
FEB_struct.Geometry.NodeSet{16}.Set = bc(4).vert_inf; 
FEB_struct.Geometry.NodeSet{16}.Name = 'vert_L4_inf';
%L5
FEB_struct.Geometry.NodeSet{17}.Set = bc(5).vert_sup; 
FEB_struct.Geometry.NodeSet{17}.Name = 'vert_L5_sup';
FEB_struct.Geometry.NodeSet{18}.Set = bc(5).vert_inf; 
FEB_struct.Geometry.NodeSet{18}.Name = 'vert_L5_inf';
%}
%Adding BC information
%fixed L2
if nbodies ==2
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;

%fixed L3
elseif nbodies ==3
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
end

%add prescribed displacement
%{
FEB_struct.Boundary.Prescribe{1}.SetName=FEB_struct.Geometry.NodeSet{5}.Name;
FEB_struct.Boundary.Prescribe{1}.bc='z';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.Scale=-1;
FEB_struct.Boundary.Prescribe{1}.Type='relative';
%}
FEB_struct.Boundary.Prescribe{1}.Set=bc(1).vert_sup;
FEB_struct.Boundary.Prescribe{1}.bc='z';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=ones(size(bc(1).vert_sup));

%Adding load information
%{
forceMagnitude=[0 0 -2000];

bcPrescribeMagnitudes=forceMagnitude(ones(1,numel(disc_sup)),:);
FEB_struct.Loads.Nodal_load{1}.bc='z';
FEB_struct.Loads.Nodal_load{1}.lc=1;
FEB_struct.Loads.Nodal_load{1}.Set=disc_sup+size(Lmes(1).VT,1)+ size(Lmes(2).VT,1)+size(Lmes(3).VT,1);
FEB_struct.Loads.Nodal_load{1}.nodeScale=bcPrescribeMagnitudes(:,1);
%}

%% CONTACTS

%contact L1-IVD1

FEB_struct.Geometry.Surface{1}.Set= F(1).vert_inf;
FEB_struct.Geometry.Surface{1}.Type='tri3';
FEB_struct.Geometry.Surface{1}.Name='Contact_master_inf_L1';
%area = area_surface(V,F(1).vert_inf);

FEB_struct.Geometry.Surface{2}.Set=FI_sup1;
FEB_struct.Geometry.Surface{2}.Type='quad4';
FEB_struct.Geometry.Surface{2}.Name='Contact_slave_sup_I1';

FEB_struct.Contact{1}.Surface{1}.SetName=FEB_struct.Geometry.Surface{1}.Name;
FEB_struct.Contact{1}.Surface{1}.Type='master';
FEB_struct.Contact{1}.Surface{2}.SetName=FEB_struct.Geometry.Surface{2}.Name;
FEB_struct.Contact{1}.Surface{2}.Type='slave';

%contact L2-IVD1
FEB_struct.Geometry.Surface{3}.Set= F(2).vert_sup;
FEB_struct.Geometry.Surface{3}.Type='tri3';
FEB_struct.Geometry.Surface{3}.Name='Contact_master_sup_L2';

FEB_struct.Geometry.Surface{4}.Set=FI_inf1;
FEB_struct.Geometry.Surface{4}.Type='quad4';
FEB_struct.Geometry.Surface{4}.Name='Contact_slave_inf_I1';

FEB_struct.Contact{2}.Surface{1}.SetName=FEB_struct.Geometry.Surface{3}.Name;
FEB_struct.Contact{2}.Surface{1}.Type='master';
FEB_struct.Contact{2}.Surface{2}.SetName=FEB_struct.Geometry.Surface{4}.Name;
FEB_struct.Contact{2}.Surface{2}.Type='slave';
if nbodies ==3

%contact L2-IVD2
    FEB_struct.Geometry.Surface{5}.Set= F(2).vert_inf;
    FEB_struct.Geometry.Surface{5}.Type='tri3';
    FEB_struct.Geometry.Surface{5}.Name='Contact_master_inf_L2';

    FEB_struct.Geometry.Surface{6}.Set=FI_sup2;
    FEB_struct.Geometry.Surface{6}.Type='quad4';
    FEB_struct.Geometry.Surface{6}.Name='Contact_slave_sup_I2';

    FEB_struct.Contact{3}.Surface{1}.SetName=FEB_struct.Geometry.Surface{5}.Name;
    FEB_struct.Contact{3}.Surface{1}.Type='master';
    FEB_struct.Contact{3}.Surface{2}.SetName=FEB_struct.Geometry.Surface{6}.Name;
    FEB_struct.Contact{3}.Surface{2}.Type='slave';

    %contact L3- IVD2
    FEB_struct.Geometry.Surface{7}.Set= F(3).vert_sup;
    FEB_struct.Geometry.Surface{7}.Type='tri3';
    FEB_struct.Geometry.Surface{7}.Name='Contact_master_sup_L3';

    FEB_struct.Geometry.Surface{8}.Set=FI_inf2;
    FEB_struct.Geometry.Surface{8}.Type='quad4';
    FEB_struct.Geometry.Surface{8}.Name='Contact_slave_inf_I2';

    FEB_struct.Contact{4}.Surface{1}.SetName=FEB_struct.Geometry.Surface{7}.Name;
    FEB_struct.Contact{4}.Surface{1}.Type='master';
    FEB_struct.Contact{4}.Surface{2}.SetName=FEB_struct.Geometry.Surface{8}.Name;
    FEB_struct.Contact{4}.Surface{2}.Type='slave';
end
%surface for pressure
%{
FEB_struct.Geometry.Surface{9}.Set= F(1).vert_sup;
FEB_struct.Geometry.Surface{9}.Type='tri3';
FEB_struct.Geometry.Surface{9}.Name='pressure';
%}


%contact L3- IVD3
%{
FEB_struct.Geometry.Surface{9}.Set= F(3).vert_inf;
FEB_struct.Geometry.Surface{9}.Type='tri3';
FEB_struct.Geometry.Surface{9}.Name='Contact_master_sup_L3';

FEB_struct.Geometry.Surface{10}.Set=FI_sup3;
FEB_struct.Geometry.Surface{10}.Type='quad4';
FEB_struct.Geometry.Surface{10}.Name='Contact_slave_sup_I3';

FEB_struct.Contact{5}.Surface{1}.SetName=FEB_struct.Geometry.Surface{9}.Name;
FEB_struct.Contact{5}.Surface{1}.Type='master';
FEB_struct.Contact{5}.Surface{2}.SetName=FEB_struct.Geometry.Surface{9}.Name;
FEB_struct.Contact{5}.Surface{2}.Type='slave';


%contact L4- IVD3
FEB_struct.Geometry.Surface{11}.Set= F(4).vert_sup;
FEB_struct.Geometry.Surface{11}.Type='tri3';
FEB_struct.Geometry.Surface{11}.Name='Contact_master_sup_L4';

FEB_struct.Geometry.Surface{12}.Set=FI_inf3;
FEB_struct.Geometry.Surface{12}.Type='quad4';
FEB_struct.Geometry.Surface{12}.Name='Contact_slave_inf_I3';

FEB_struct.Contact{6}.Surface{1}.SetName=FEB_struct.Geometry.Surface{11}.Name;
FEB_struct.Contact{6}.Surface{1}.Type='master';
FEB_struct.Contact{6}.Surface{2}.SetName=FEB_struct.Geometry.Surface{12}.Name;
FEB_struct.Contact{6}.Surface{2}.Type='slave';
%contact L4- IVD4
FEB_struct.Geometry.Surface{13}.Set= F(4).vert_inf;
FEB_struct.Geometry.Surface{13}.Type='tri3';
FEB_struct.Geometry.Surface{13}.Name='Contact_master_inf_L4';

FEB_struct.Geometry.Surface{14}.Set=FI_sup4;
FEB_struct.Geometry.Surface{14}.Type='quad4';
FEB_struct.Geometry.Surface{14}.Name='Contact_slave_sup_I4';

FEB_struct.Contact{7}.Surface{1}.SetName=FEB_struct.Geometry.Surface{13}.Name;
FEB_struct.Contact{7}.Surface{1}.Type='master';
FEB_struct.Contact{7}.Surface{2}.SetName=FEB_struct.Geometry.Surface{14}.Name;
FEB_struct.Contact{7}.Surface{2}.Type='slave';
%contact L5- IVD4
FEB_struct.Geometry.Surface{15}.Set= F(5).vert_sup;
FEB_struct.Geometry.Surface{15}.Type='tri3';
FEB_struct.Geometry.Surface{15}.Name='Contact_master_sup_L5';

FEB_struct.Geometry.Surface{16}.Set=FI_inf4;
FEB_struct.Geometry.Surface{16}.Type='quad4';
FEB_struct.Geometry.Surface{16}.Name='Contact_slave_inf_I5';

FEB_struct.Contact{8}.Surface{1}.SetName=FEB_struct.Geometry.Surface{15}.Name;
FEB_struct.Contact{8}.Surface{1}.Type='master';
FEB_struct.Contact{8}.Surface{2}.SetName=FEB_struct.Geometry.Surface{16}.Name;
FEB_struct.Contact{8}.Surface{2}.Type='slave';
%}

%% add discreteset for ligaments



%% contact definition
FEB_struct.Contact{1}.Type='tied';
FEB_struct.Contact{1}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{1}.Values={0,tolerance,contactPenalty,0,10,0.01};

FEB_struct.Contact{2}.Type='tied';
FEB_struct.Contact{2}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{2}.Values={0,tolerance,contactPenalty,0,10,0.01};
if nbodies==3
    FEB_struct.Contact{3}.Type='tied';
    FEB_struct.Contact{3}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
    FEB_struct.Contact{3}.Values={0,tolerance,contactPenalty,0,10,0.01};

    FEB_struct.Contact{4}.Type='tied';
    FEB_struct.Contact{4}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
    FEB_struct.Contact{4}.Values={0,tolerance,contactPenalty,0,10,0.01};
end
%{
FEB_struct.Contact{5}.Type='tied';
FEB_struct.Contact{5}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{5}.Values={0,tolerance,contactPenalty,0,10,0.01};

FEB_struct.Contact{6}.Type='tied';
FEB_struct.Contact{6}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{6}.Values={0,tolerance,contactPenalty,0,10,0.01};

FEB_struct.Contact{7}.Type='tied';
FEB_struct.Contact{7}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{7}.Values={0,tolerance,contactPenalty,0,10,0.01};

FEB_struct.Contact{8}.Type='tied';
FEB_struct.Contact{8}.Properties={'laugon','tolerance','penalty','minaug','maxaug','search_tolerance'};
FEB_struct.Contact{8}.Values={0,tolerance,contactPenalty,0,10,0.01};
%}

%% Control section
FEB_struct.Step{1}.Control.AnalysisType='static';
FEB_struct.Step{1}.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Step{1}.Control.Values={10,0.1,15,0,0.001,0.01,0,0.9};
FEB_struct.Step{1}.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter'};
FEB_struct.Step{1}.Control.TimeStepperValues={1e-5,0.1,10,10};
%FEB_struct.Step{2}.Control=FEB_struct.Step{1}.Control;


%Specify log file output
run_output_name1=[modelName,'_node_out.txt'];
run_output_name2=[modelName,'_F_out.txt'];
run_output_name3=[modelName,'_R_out.txt'];
run_output_name4=[modelName,'_S_out.txt'];
FEB_struct.run_output_names={run_output_name1,run_output_name2,run_output_name3,run_output_name4};
FEB_struct.output_types={'node_data','element_data','node_data','element_data',};
FEB_struct.data_types={'ux;uy;uz','Fxx;Fxy;Fxz;Fyx;Fyy;Fyz;Fzx;Fzy;Fzz','Rx;Ry;Rz','sx;sy;sz'};
%fprintf(fid,'<element_data data="Ex;Ey;Ez" name="element strain"file="bdyn_%d_node_strain_out_linear_iso.txt"> </element_data>\n',j);
%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume'};

%Load curves
FEB_struct.LoadData.LoadCurves.id=[1 2];
FEB_struct.LoadData.LoadCurves.type={'linear' 'linear' 'linear' 'linear' ...
'linear' 'linear' 'linear' 'linear' 'linear' };
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 -1];[0 0;1 0.2]};

%{
%pressure 2000N over "area"
area = area_surface(V,F(1).vert_sup)
pressureLevel = 2000/area %N/mm^2

%Adding load information

FEB_struct.Loads.Surface_load{1}.Type='pressure';
% FEB_struct.Loads.Surface_load{1}.Set=F_pressure;
FEB_struct.Loads.Surface_load{1}.SetName=FEB_struct.Geometry.Surface{9}.Name;
FEB_struct.Loads.Surface_load{1}.lcPar='pressure';
FEB_struct.Loads.Surface_load{1}.lcParValue=pressureLevel;
FEB_struct.Loads.Surface_load{1}.lc=1;
%}
%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Turn on displaying of progress

febStruct2febFile(FEB_struct);

%% run FEBio
 %% run FEBio
    FEBioRunStruct.run_filename=FEB_struct.run_filename;
    FEBioRunStruct.run_logname=FEB_struct.run_logname;
    FEBioRunStruct.disp_on=1;
    FEBioRunStruct.disp_log_on=1;
    FEBioRunStruct.runMode= 'internal';%'external';
    FEBioRunStruct.t_check=5; %Time for checking log file (dont set too small)
    FEBioRunStruct.maxtpi=1e99; %Max analysis time
    FEBioRunStruct.maxLogCheckTime=10; %Max log file checking time
    %FEBioRunStruct.FEBioPath='C:\Program Files\febio-2.5.1\bin\febio2.exe';
    
   [runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!
      %{
    if runFlag==1
        %savepath = sprintf('C:\Users\lavecchc\Desktop\sensitivity');
 %       fullfile(savePath,FEB_struct.run_output_names{1})
    [~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements

    DN=N_disp_mat(FI_sup1,end,:); %Final nodal displacements
    
    %max_disp(j,:) = min(DN); %because the displacement is negative
    
    V_def=VT+DN;
    DN_magnitude=sqrt(sum(DN.^2,2));
    
    define fb
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);

    hf1=cFigure;
    title('The deformed model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

    hps=patch('Faces',Fbnew,'Vertices',VT,'FaceColor','flat','CData',CF);

    view(3); axis tight;  axis equal;  grid on;
    colormap gjet; colorbar;
    % camlight headlight;
    set(gca,'FontSize',fontSize);
    drawnow;
   
    [time_mat, N_force_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{2}); %Nodal forces
    time_mat=[0; time_mat(:)]; %Time
    
    %DERIVING STRESS METRICS

    %Get Z forces
    FZ=sum(N_force_mat(FI_sup1,end,:),1);
    FZ=[0; FZ(:)]; %Mean top surface nodal forces

    %Derive applied stretch
    DZ_set=N_disp_mat(FI_sup1,end,:); %Final nodal displacements
    DZ_set=mean(DZ_set,1);
    %ezample
    sampleHeight = 10;
    
    stretch_sim=(DZ_set+sampleHeight)./sampleHeight;
    stretch_sim=[1; stretch_sim(:)];
    
    initialArea = 20;
    %Derive simulated Cauchy stress (alternatively import stress and take the mean)
    currentArea=initialArea./stretch_sim;
    stress_cauchy_sim=FZ./currentArea; %Cauchy stress
    stress_cauchy_sim=stress_cauchy_sim.*1e3; %Scale to kPa
    hf1=cFigure;
    title('Stress curve, upper plate L1','FontSize',fontSize);
    xlabel('Time (s)','FontSize',fontSize); ylabel('\sigma Cauchy stress (kPa)','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    lineWidth = 1;
    plot(time_mat(:),stress_cauchy_sim(:),'r.-','lineWidth',lineWidth,'markerSize',markerSize);

    view(2); axis tight;  grid on;
    set(gca,'FontSize',fontSize);
    drawnow;
    
    end

%}

 end
    %}