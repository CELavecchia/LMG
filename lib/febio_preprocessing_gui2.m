%% FEBio preprocessing - FU L1-L2-L3
function FEB_struct = febio_preprocessing_gui(modelName, nbodies, Lmes, mesh_struct_IVD2, runFebio)

%%
% fibres embedded

%clear; close all; clc;
%% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath));



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


% contact properties
contactPenalty=1;
tolerance =  0.5;
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


%% import datasets

%==============================================================
fprintf('import vertebrae dataset\n\n');

% IVD 1-2
VT = mesh_struct_IVD2(1).V;
E_Ivd1 = mesh_struct_IVD2(1).E;
bC = mesh_struct_IVD2(1).faceBoundaryMarker;
Fbnew = mesh_struct_IVD2(1).Fb;
C = mesh_struct_IVD2(1).C;
F1 = mesh_struct_IVD2(1).FE;

% IVD 2-3
VT2 = mesh_struct_IVD2(2).V;
E_Ivd2 = mesh_struct_IVD2(2).E;
bC2 = mesh_struct_IVD2(2).faceBoundaryMarker;
Fbnew2 = mesh_struct_IVD2(2).Fb;
C2 = mesh_struct_IVD2(2).C;
F2 = mesh_struct_IVD2(2).FE;

% IVD 3-4
VT3 = mesh_struct_IVD2(3).V;
E_Ivd3 = mesh_struct_IVD2(3).E;
bC3 = mesh_struct_IVD2(3).faceBoundaryMarker;
Fbnew3 = mesh_struct_IVD2(3).Fb;
C3 = mesh_struct_IVD2(3).C;
F3 = mesh_struct_IVD2(3).FE;

% IVD 4-5
VT4 = mesh_struct_IVD2(4).V;
E_Ivd4 = mesh_struct_IVD2(4).E;
bC4 = mesh_struct_IVD2(4).faceBoundaryMarker;
Fbnew4 = mesh_struct_IVD2(4).Fb;
C4 = mesh_struct_IVD2(4).C;
F4 = mesh_struct_IVD2(4).FE;



%% define nodeset for the attachment points of the ligaments  ----- move this function!!
%ligam = attachment_points(Lmes);


%% join datasets

%[x_node_IVD1, y_node_IVD1, x_node_IVD2, y_node_IVD2] = nodeset(VT,VT2);

if nbodies ==2              % FU L1 - IVD1 - L2
    
    %--------------- only to avoid the initial penetration-----------
    %Lmes(2).VT(:,3) = Lmes(2).VT(:,3)-5;
    %VT(:,3) = VT(:,3) -2.5;
    
    cFigure;
    
    subplot(1,2,1)
    title('FEBio model','FontSize',fontSize);
    patch('Faces',Lmes(1).Fb,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Lmes(1).Cb,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    patch('Faces',Lmes(2).Fb,'Vertices',Lmes(2).VT,'FaceColor','flat','CData',Lmes(2).Cb,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    patch('Faces',F1,'Vertices',VT,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    axis tight;  axis equal;  grid on;
    colormap(autumn);
    camlight headlight;
    set(gca,'FontSize',fontSize);
    
    subplot(1,2,2)
    title('FEBio pre-processing-BC applied','FontSize',fontSize);
    patch('Faces',Lmes(1).Fb,'Vertices',Lmes(1).VT,'FaceColor','flat','CData',Lmes(1).Cb,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    patch('Faces',Lmes(2).Fb,'Vertices',Lmes(2).VT,'FaceColor','flat','CData',Lmes(2).Cb,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    patch('Faces',F1,'Vertices',VT,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    axis tight;  axis equal;  grid on;hold on;
    
    %-------------------------------------------------------
    
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
    
elseif nbodies ==5          %  L1 - IVD1 - L2 - IVD2 - L3 - IVD3 - L4 - IVD4 - L5
    V = [Lmes(1).VT; Lmes(2).VT; Lmes(3).VT; Lmes(4).VT;Lmes(5).VT; VT; VT2; VT3; VT4];
    E1 = Lmes(1).E;
    E2 = Lmes(2).E + size(Lmes(1).VT,1);
    E3 = Lmes(3).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1);
    E4 = Lmes(4).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1) + size(Lmes(3).VT,1);
    E5 = Lmes(5).E + size(Lmes(1).VT,1)+ size(Lmes(2).VT,1) + size(Lmes(3).VT,1)+ size(Lmes(4).VT,1);
    Lmes(2).Fb = Lmes(2).Fb + size(Lmes(1).VT,1);
    Lmes(3).Fb = Lmes(3).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1);
    Lmes(4).Fb = Lmes(4).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1)+size(Lmes(3).VT,1);
    Lmes(5).Fb = Lmes(5).Fb + size(Lmes(1).VT,1)+size(Lmes(2).VT,1)+size(Lmes(3).VT,1)+size(Lmes(4).VT,1);
    Fbnew = Fbnew + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) ;%IVD1
    Fbnew2 = Fbnew2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1);%IVD2
    Fbnew3 = Fbnew3 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1) + size(VT,1)+ size(VT2,1);%IVD3
    Fbnew4 = Fbnew4 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+...
        size(Lmes(4).VT,1) + size(VT,1)+ size(VT2,1)+ size(VT3,1);%IVD3
    Fbnew4 = Fbnew4 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ ...
        size(Lmes(4).VT,1) + size(Lmes(5).VT,1) + size(VT,1)+ size(VT2,1)+ size(VT3,1)+ size(VT4,1);%IVD3
    
    
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
    
    febMatL5= Lmes(5).elementMaterialIndices;
    febMatL5(Lmes(5).elementMaterialIndices==-2)=1;
    febMatL5(Lmes(5).elementMaterialIndices==-3)=2;
    
    E_ivd1 = E_Ivd1 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1);
    E_ivd1 = [E_ivd1(:,5:8) E_ivd1(:,1:4)]; %it was inverted for febio
    E_ivd2 = E_Ivd2 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1)+ size(VT,1)+ size(VT2,1);
    E_ivd2 = [E_ivd2(:,5:8) E_ivd2(:,1:4)];%it was inverted for febio
    E_ivd3 = E_Ivd3 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ size(Lmes(4).VT,1)+ size(VT,1)+ size(VT2,1)+ size(VT3,1);
    E_ivd3 = [E_ivd3(:,5:8) E_ivd3(:,1:4)];%it was inverted for febio
    E_ivd4 = E_Ivd4 + size(Lmes(1).VT,1) + size(Lmes(2).VT,1)+ size(Lmes(3).VT,1)+ ...
        size(Lmes(4).VT,1)+ size(Lmes(5).VT,1)+ size(VT,1)+ size(VT2,1)+ size(VT3,1)+ size(VT4,1);
    E_ivd4 = [E_ivd4(:,5:8) E_ivd4(:,1:4)];%it was inverted for febio
    
    
    
end

%nodeset for the output
%{
x_node_IVD1 =  x_node_IVD1 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1);
y_node_IVD1 = y_node_IVD1 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1);
x_node_IVD2 = x_node_IVD2 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1) +size(VT(:,1),1);
y_node_IVD2 = y_node_IVD2 + size(Lmes(1).VT(:,1),1) +size(Lmes(2).VT(:,1),1)+ size(Lmes(3).VT(:,1),1) +size(VT(:,1),1);
%}


% update indices of the ligaments attaching points
%ligam2 = update_indeces(ligam, size(Lmes(1).VT(:,1),1), size(Lmes(2).VT(:,1)) );


%element indices for the IVD
febMatIVD= C;
%anulus   % you can select also different material prop for each layer
febMatIVD(C==2)=3;%3; %nucleus
febMatIVD(C==3)=3;%4;
febMatIVD(C==4)=3;%3;
febMatIVD(C==5)=3;%4;
febMatIVD(C==6)=3;%3;
febMatIVD(C==7)=3;%4;
febMatIVD(C==8)=3;%3;
febMatIVD(C==9)=3;%4;
%nucleus
febMatIVD(C==10)=3;%5;
febMatIVD(C==11)=3;%5;


%IVD2
febMatIVD2= C2;
febMatIVD2(C2==2)=3;%3; %nucleus
febMatIVD2(C2==3)=3;%4;
febMatIVD2(C2==4)=3;%3;
febMatIVD2(C2==5)=3;%4;
febMatIVD2(C2==6)=3;%3;
febMatIVD2(C2==7)=3;%4;
febMatIVD2(C2==8)=3;%3;
febMatIVD2(C2==9)=3;%4;
%nucleus
febMatIVD2(C2==10)=3;%5;
febMatIVD2(C2==11)=3;%5;

% IVD3
febMatIVD3= C3;
febMatIVD3(C3==2)=3;%3; %nucleus
febMatIVD3(C3==3)=3;%4;
febMatIVD3(C3==4)=3;%3;
febMatIVD3(C3==5)=3;%4;
febMatIVD3(C3==6)=3;%3;
febMatIVD3(C3==7)=3;%4;
febMatIVD3(C3==8)=3;%3;
febMatIVD3(C3==9)=3;%4;
%nucleus
febMatIVD3(C3==10)=3;%5;
febMatIVD3(C3==11)=3;%5;

% IVD4
febMatIVD4= C4;
febMatIVD4(C4==2)=3;%3; %nucleus
febMatIVD4(C4==3)=3;%4;
febMatIVD4(C4==4)=3;%3;
febMatIVD4(C4==5)=3;%4;
febMatIVD4(C4==6)=3;%3;
febMatIVD4(C4==7)=3;%4;
febMatIVD4(C4==8)=3;%3;
febMatIVD4(C4==9)=3;%4;
%nucleus
febMatIVD4(C4==10)=3;%5;
febMatIVD4(C4==11)=3;%5;


%% CREATING FIBRE DIRECTIONS
% check other script
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
    model.E3 = E3;model.E_ivd2 = E_ivd2; %elements
    model.matind3 =febMatL3; %mat ind vert
    model.FbI2 = febMatIVD2;    %mat ind IVD
    model.bcI2 = mesh_struct_IVD2(2).faceBoundaryMarker; %bc IVD
    model.Fb3 = Lmes(3).Fb; %faces vert
    model.Cb3 = Lmes(3).Cb; % cb vert
    model.E4 = E4;model.E_ivd3 = E_ivd3; %elements
    model.matind4 =febMatL4; %mat ind vert
    model.FbI3 = febMatIVD3;    %mat ind IVD
    model.bcI3 = mesh_struct_IVD2(3).faceBoundaryMarker; %bc IVD
    model.Fb4 = Lmes(4).Fb; %faces vert
    model.Cb4 = Lmes(4).Cb; % cb vert
elseif nbodies ==5
    model.E3 = E3;model.E_ivd2 = E_ivd2; %elements
    model.matind3 =febMatL3; %mat ind vert
    model.FbI2 = febMatIVD2;    %mat ind IVD
    model.bcI2 = mesh_struct_IVD2(2).faceBoundaryMarker; %bc IVD
    model.Fb3 = Lmes(3).Fb; %faces vert
    model.Cb3 = Lmes(3).Cb; % cb vert
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
bcIVD(1).sup=unique(FI_sup1(:)); %disc_sup1
logicFace=bC==1;
FI_inf1=Fbnew(logicFace,:);
bcIVD(1).inf=unique(FI_inf1(:)); %disc_inf1

logicFace=bC2==2;
FI_sup2=Fbnew2(logicFace,:);
bcIVD(2).sup=unique(FI_sup2(:));
logicFace=bC2==1;
FI_inf2=Fbnew2(logicFace,:);
bcIVD(2).inf=unique(FI_inf2(:));

logicFace=bC3==1;
FI_sup3=Fbnew3(logicFace,:);
bcIVD(3).sup=unique(FI_sup3(:));
logicFace=bC3==2;
FI_inf3=Fbnew3(logicFace,:);
bcIVD(3).inf=unique(FI_inf2(:));

logicFace=bC4==1;
FI_sup4=Fbnew4(logicFace,:);
bcIVD(4).sup=unique(FI_sup4(:));
logicFace=bC4==2;
FI_inf4=Fbnew4(logicFace,:);
bcIVD(4).inf=unique(FI_inf4(:));


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
elseif nbodies==4
    FEB_struct.Geometry.Elements={E1 E2 E3 E4 E_ivd1 E_ivd2 E_ivd3}; %The element sets
    FEB_struct.Geometry.ElementType={'tet4','tet4','tet4','tet4','hex8','hex8','hex8'}; %The element types
    FEB_struct.Geometry.ElementMat={febMatL1', febMatL2',febMatL3',febMatL4', febMatIVD', febMatIVD2',febMatIVD3'};
    FEB_struct.Geometry.ElementsPartName={'L1' 'L2' 'L3' 'L4' 'IVD1' 'IVD2' 'IVD3'};
elseif nbodies==5
    FEB_struct.Geometry.Elements={E1 E2 E3 E4 E5 E_ivd1 E_ivd2 E_ivd3 E_ivd4}; %The element sets
    FEB_struct.Geometry.ElementType={'tet4','tet4','tet4','tet4','tet4','hex8','hex8','hex8','hex8'}; %The element types
    FEB_struct.Geometry.ElementMat={febMatL1', febMatL2',febMatL3',febMatL4',febMatL5', febMatIVD', febMatIVD2',febMatIVD3',febMatIVD4'};
    FEB_struct.Geometry.ElementsPartName={'L1' 'L2' 'L3' 'L4' 'L5' 'IVD1' 'IVD2' 'IVD3' 'IVD4'};
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




%IVD   ---- only for this test, check in the other script for the material
%properties for AF, NP and CEP
%isotropic linear
FEB_struct.Materials{3}.Type = 'isotropic elastic';
FEB_struct.Materials{3}.Name = 'nucleus'; %np1
FEB_struct.Materials{3}.Properties = {'E','v'};
FEB_struct.Materials{3}.Values = {Eanul,vanul};



%CEP
%check this one... the elements have been selected
%{
FEB_struct.Materials{6}.Type = 'isotropic elastic';
FEB_struct.Materials{6}.Name = 'bone';
FEB_struct.Materials{6}.Properties = {'E','v'};
FEB_struct.Materials{6}.Values = {Ecep,vcep};
%}
%% FIBRES DEFINITION ---> RE-CHECK!



%% BCs

%Defining node sets


%IVD1
FEB_struct.Geometry.NodeSet{1}.Set = bcIVD(1).sup;
FEB_struct.Geometry.NodeSet{1}.Name = 'disc_sup1';
FEB_struct.Geometry.NodeSet{2}.Set = bcIVD(1).inf;
FEB_struct.Geometry.NodeSet{2}.Name = 'disc_inf1';

if nbodies==2
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
    
    % BC
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{6}.Name;
    
    
    plotV(V(bc(2).vert_inf,:),'k.','MarkerSize',markerSize),hold on;
    %  colormap(autumn);
    camlight headlight;
    set(gca,'FontSize',fontSize);
    
    
elseif nbodies==3
    %IVD2
    FEB_struct.Geometry.NodeSet{3}.Set = bcIVD(2).sup;
    FEB_struct.Geometry.NodeSet{3}.Name = 'disc_sup2';
    FEB_struct.Geometry.NodeSet{4}.Set = bcIVD(2).sup;
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
    
    % bc
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{10}.Name;
    
elseif nbodies==4
    %IVD2
    FEB_struct.Geometry.NodeSet{3}.Set = bcIVD(2).sup;
    FEB_struct.Geometry.NodeSet{3}.Name = 'disc_sup2';
    FEB_struct.Geometry.NodeSet{4}.Set = bcIVD(2).sup;
    FEB_struct.Geometry.NodeSet{4}.Name = 'disc_inf2';
    %IVD3
    FEB_struct.Geometry.NodeSet{5}.Set = bcIVD(3).sup;
    FEB_struct.Geometry.NodeSet{5}.Name = 'disc_sup3';
    FEB_struct.Geometry.NodeSet{6}.Set = bcIVD(3).sup;
    FEB_struct.Geometry.NodeSet{6}.Name = 'disc_inf3';
    %L1
    FEB_struct.Geometry.NodeSet{7}.Set = bc(1).vert_sup;
    FEB_struct.Geometry.NodeSet{7}.Name = 'vert_L1_sup';
    FEB_struct.Geometry.NodeSet{8}.Set = bc(1).vert_inf;
    FEB_struct.Geometry.NodeSet{8}.Name = 'vert_L1_inf';
    %L2
    FEB_struct.Geometry.NodeSet{9}.Set = bc(2).vert_sup;
    FEB_struct.Geometry.NodeSet{9}.Name = 'vert_L2_sup';
    FEB_struct.Geometry.NodeSet{10}.Set = bc(2).vert_inf;
    FEB_struct.Geometry.NodeSet{10}.Name = 'vert_L2_inf';
    %L3
    FEB_struct.Geometry.NodeSet{11}.Set = bc(3).vert_sup;
    FEB_struct.Geometry.NodeSet{11}.Name = 'vert_L3_sup';
    FEB_struct.Geometry.NodeSet{12}.Set = bc(3).vert_inf;
    FEB_struct.Geometry.NodeSet{12}.Name = 'vert_L3_inf';
    %L4
    FEB_struct.Geometry.NodeSet{13}.Set = bc(4).vert_sup;
    FEB_struct.Geometry.NodeSet{13}.Name = 'vert_L4_sup';
    FEB_struct.Geometry.NodeSet{14}.Set = bc(4).vert_inf;
    FEB_struct.Geometry.NodeSet{14}.Name = 'vert_L4_inf';
    
    % BC
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{14}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{14}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{14}.Name;
    
    
elseif nbodies==5
    %IVD2
    FEB_struct.Geometry.NodeSet{3}.Set = bcIVD(2).sup;
    FEB_struct.Geometry.NodeSet{3}.Name = 'disc_sup2';
    FEB_struct.Geometry.NodeSet{4}.Set = bcIVD(2).sup;
    FEB_struct.Geometry.NodeSet{4}.Name = 'disc_inf2';
    %IVD3
    FEB_struct.Geometry.NodeSet{5}.Set = bcIVD(3).sup;
    FEB_struct.Geometry.NodeSet{5}.Name = 'disc_sup3';
    FEB_struct.Geometry.NodeSet{6}.Set = bcIVD(3).sup;
    FEB_struct.Geometry.NodeSet{6}.Name = 'disc_inf3';
    %IVD4
    FEB_struct.Geometry.NodeSet{7}.Set = bcIVD(4).sup;
    FEB_struct.Geometry.NodeSet{7}.Name = 'disc_sup4';
    FEB_struct.Geometry.NodeSet{8}.Set = bcIVD(4).sup;
    FEB_struct.Geometry.NodeSet{8}.Name = 'disc_inf4';
    %L1
    FEB_struct.Geometry.NodeSet{9}.Set = bc(1).vert_sup;
    FEB_struct.Geometry.NodeSet{9}.Name = 'vert_L1_sup';
    FEB_struct.Geometry.NodeSet{10}.Set = bc(1).vert_inf;
    FEB_struct.Geometry.NodeSet{10}.Name = 'vert_L1_inf';
    %L2
    FEB_struct.Geometry.NodeSet{11}.Set = bc(2).vert_sup;
    FEB_struct.Geometry.NodeSet{11}.Name = 'vert_L2_sup';
    FEB_struct.Geometry.NodeSet{12}.Set = bc(2).vert_inf;
    FEB_struct.Geometry.NodeSet{12}.Name = 'vert_L2_inf';
    %L3
    FEB_struct.Geometry.NodeSet{13}.Set = bc(3).vert_sup;
    FEB_struct.Geometry.NodeSet{13}.Name = 'vert_L3_sup';
    FEB_struct.Geometry.NodeSet{14}.Set = bc(3).vert_inf;
    FEB_struct.Geometry.NodeSet{14}.Name = 'vert_L3_inf';
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
    
    % BC
    FEB_struct.Boundary.Fix{1}.bc='x';
    FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{18}.Name;
    FEB_struct.Boundary.Fix{2}.bc='y';
    FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{18}.Name;
    FEB_struct.Boundary.Fix{3}.bc='z';
    FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{18}.Name;
    
end


%add prescribed displacement

FEB_struct.Boundary.Prescribe{1}.Set=bc(1).vert_sup;
FEB_struct.Boundary.Prescribe{1}.bc='z';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=ones(size(bc(1).vert_sup));

plotV(V(bc(1).vert_sup,:),'b.','MarkerSize',markerSize),hold on;
%colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
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
    
elseif  nbodies ==4
    
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
    %contact L4- IVD3
    FEB_struct.Geometry.Surface{9}.Set= F(4).vert_sup;
    FEB_struct.Geometry.Surface{9}.Type='tri3';
    FEB_struct.Geometry.Surface{9}.Name='Contact_master_sup_L4';
    
    FEB_struct.Geometry.Surface{10}.Set=FI_inf3;
    FEB_struct.Geometry.Surface{10}.Type='quad4';
    FEB_struct.Geometry.Surface{10}.Name='Contact_slave_inf_I3';
    
    FEB_struct.Contact{5}.Surface{1}.SetName=FEB_struct.Geometry.Surface{9}.Name;
    FEB_struct.Contact{5}.Surface{1}.Type='master';
    FEB_struct.Contact{5}.Surface{2}.SetName=FEB_struct.Geometry.Surface{10}.Name;
    FEB_struct.Contact{5}.Surface{2}.Type='slave';
    
elseif  nbodies ==5
    
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
    %contact L4- IVD3
    FEB_struct.Geometry.Surface{9}.Set= F(4).vert_sup;
    FEB_struct.Geometry.Surface{9}.Type='tri3';
    FEB_struct.Geometry.Surface{9}.Name='Contact_master_sup_L4';
    
    FEB_struct.Geometry.Surface{10}.Set=FI_inf3;
    FEB_struct.Geometry.Surface{10}.Type='quad4';
    FEB_struct.Geometry.Surface{10}.Name='Contact_slave_inf_I3';
    
    FEB_struct.Contact{5}.Surface{1}.SetName=FEB_struct.Geometry.Surface{9}.Name;
    FEB_struct.Contact{5}.Surface{1}.Type='master';
    FEB_struct.Contact{5}.Surface{2}.SetName=FEB_struct.Geometry.Surface{10}.Name;
    FEB_struct.Contact{5}.Surface{2}.Type='slave';
    %contact L5- IVD4
    FEB_struct.Geometry.Surface{11}.Set= F(5).vert_sup;
    FEB_struct.Geometry.Surface{11}.Type='tri3';
    FEB_struct.Geometry.Surface{11}.Name='Contact_master_sup_L5';
    
    FEB_struct.Geometry.Surface{12}.Set=FI_inf4;
    FEB_struct.Geometry.Surface{12}.Type='quad4';
    FEB_struct.Geometry.Surface{12}.Name='Contact_slave_inf_I4';
    
    FEB_struct.Contact{6}.Surface{1}.SetName=FEB_struct.Geometry.Surface{11}.Name;
    FEB_struct.Contact{6}.Surface{1}.Type='master';
    FEB_struct.Contact{6}.Surface{2}.SetName=FEB_struct.Geometry.Surface{12}.Name;
    FEB_struct.Contact{6}.Surface{2}.Type='slave';
end

%surface for pressure
%{
FEB_struct.Geometry.Surface{9}.Set= F(1).vert_sup;
FEB_struct.Geometry.Surface{9}.Type='tri3';
FEB_struct.Geometry.Surface{9}.Name='pressure';
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
    FEB_struct.Contact{3} = FEB_struct.Contact{1};
    FEB_struct.Contact{4} = FEB_struct.Contact{1};    
elseif nbodies==4
    FEB_struct.Contact{3} = FEB_struct.Contact{1};
    FEB_struct.Contact{4} = FEB_struct.Contact{1};
    FEB_struct.Contact{5} = FEB_struct.Contact{1};
    FEB_struct.Contact{6} = FEB_struct.Contact{1};
elseif nbodies==5
    FEB_struct.Contact{3} = FEB_struct.Contact{1};
    FEB_struct.Contact{4} = FEB_struct.Contact{1};
    FEB_struct.Contact{5} = FEB_struct.Contact{1};
    FEB_struct.Contact{6} = FEB_struct.Contact{1};
    FEB_struct.Contact{7} = FEB_struct.Contact{6};
    FEB_struct.Contact{8} = FEB_struct.Contact{6};    
end

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
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1;
FEBioRunStruct.disp_log_on=1;
FEBioRunStruct.runMode= 'internal';%'external';
FEBioRunStruct.t_check=5; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=10; %Max log file checking time
%FEBioRunStruct.FEBioPath='C:\Program Files\febio-2.5.1\bin\febio2.exe';


if runFebio == 1 
    [runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!
end
    
    %{
if runFlag==1
      
 %% POST-PROCESSING
    [time_mat, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements
    
    DN=N_disp_mat(:,2:end,:);
    DN_magnitude=sqrt(sum(DN.^2,2));
    %nodeset in deformed state
    V_def=V+DN;
    
    %Visualise the results
    hf1=cFigure;
    title('The deformed model (Z)','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

    patch('Faces',Fbnew,'Vertices',V_def(:,:,end),'FaceColor','flat','CData',DN_magnitude(:,:,end));
    patch('Faces',Fbnew2,'Vertices',V_def(:,:,end),'FaceColor','flat','CData',DN_magnitude(:,:,end));
    patch('Faces',Lmes(1).Fb,'Vertices',V_def(:,:,end),'FaceColor','flat','CData',DN_magnitude(:,:,end));
    patch('Faces',Lmes(2).Fb,'Vertices',V_def(:,:,end),'FaceColor','flat','CData',DN_magnitude(:,:,end));
    patch('Faces',Lmes(3).Fb,'Vertices',V_def(:,:,end),'FaceColor','flat','CData',DN_magnitude(:,:,end));

    view; grid on;
    colormap jet; colorbar;
%}
end
