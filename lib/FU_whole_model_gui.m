    % working on the FE model of L1-L2 and IVD

%clear; close all; clc;
function [Lmes] = FU_whole_model_gui(stl)
%%
nbodies =5; %mesh all the vertebrae
% Plot settings
figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=50; 

% path names
filePath=mfilename('fullpath');

%path to get the datasets
defaultFolder = fileparts(mfilename('fullpath'));
inputpathName=fullfile(defaultFolder,'data','output','point_cloud'); 

%path to save the dataset
outputpathName=fullfile(defaultFolder,'data','output','mat'); 
outputName=fullfile(outputpathName,'Lmes_v2.mat');

%%
%==============================================================
fprintf('import vertebrae point clouds to mesh \n');

%-----------------      import point clouds          ----------------------
for(j = 1:5)
    
    disp(sprintf('\n\n Working on the vertenbra L%d \n',j))
    
    name_body=fullfile(inputpathName,sprintf('L_body_%d.txt',j));
    %name_body = sprintf('./output/point_cloud/L_body_%d.txt',j);
    L_body = load(name_body);
    
    name_lam = fullfile(inputpathName,sprintf('L_lam_%d.txt',j));%sprintf('./output/point_cloud/L_lam_%d.txt',j);
    L_lam = load(name_lam);
    
    name_ped = fullfile(inputpathName,sprintf('L_ped_%d.txt',j));%sprintf('./output/point_cloud/L_ped_%d.txt',j);
    L_ped = load(name_ped);
    
    name_proc = fullfile(inputpathName,sprintf('L_proc_grid_%d.txt',j));%sprintf('./output/point_cloud/L_proc_grid_%d.txt',j);
    L_proc = load(name_proc);

%----------------------------------------------------------------------
%divide in two different dataasedts, left (sx) and rights(dx) ones
[L_proc_dx,L_proc_sx, L_ped_dx, L_ped_sx] = divide_dx_sx_mesh(L_proc, L_ped);

%------ remove the overlappingpoints 

%vertebral body
fprintf('------------- meshing the v body  ------------\n');
[body, Eb,Vb] = obtain_geometry_gui( L_body);
Eb = body.Mesh.Elements';
Vb = body.Mesh.Nodes';
[Fb,Cb] = element2patch(Eb,'tet4');
%figure
%patch('Faces',Fb,'Vertices',Vb,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');

fprintf('------------- meshing the pedicles and processes  ------------\n');

%processes
[ proc_ped_dx,Epp_dx,Vpp_dx] = obtain_geometry2_gui( L_proc_dx(:,1:3));
Epp_dx = proc_ped_dx.Mesh.Elements';
Vpp_dx = proc_ped_dx.Mesh.Nodes';
[Fpp_dx,Cpp_dx] = element2patch(Epp_dx,'tet4');
%figure
%patch('Faces',Fpp_dx,'Vertices',Vpp_dx,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');


[ proc_ped_sx, Epp_sx,Vpp_sx] = obtain_geometry2_gui( L_proc_sx(:,1:3));
Epp_sx = proc_ped_sx.Mesh.Elements';
Vpp_sx = proc_ped_sx.Mesh.Nodes';
[Fpp_sx,Cpp_sx] = element2patch(Epp_sx,'tet4');
%figure
%patch('Faces',Fpp_sx,'Vertices',Vpp_sx,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');


%pedicles
[ ped_dx,Epe_dx,Vpe_dx] = obtain_geometry3_gui(  L_ped_dx(:,1:3));
Epe_dx = ped_dx.Mesh.Elements';
Vpe_dx = ped_dx.Mesh.Nodes';
[Fpe_dx,Cpe_dx] = element2patch(Epe_dx,'tet4');
%figure
%patch('Faces',Fpe_dx,'Vertices',Vpe_dx,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');


[ ped_sx, Epe_sx,Vpe_sx] = obtain_geometry3_gui(  L_ped_sx(:,1:3));
Epe_sx = ped_sx.Mesh.Elements';
Vpe_sx = ped_sx.Mesh.Nodes';   
[Fpe_sx,Cpe_sx] = element2patch(Epe_sx,'tet4');
%figure
%patch('Faces',Fpe_sx,'Vertices',Vpe_sx,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');

%lamina

fprintf('------------- meshing the lamina  ------------\n');
%figure
[lam, El,Vl] = obtain_geometry_lam_gui(L_lam);
El = lam.Mesh.Elements';
Vl = lam.Mesh.Nodes';
[Fl,Cl] = element2patch(El,'tet4');
%patch('Faces',Fl,'Vertices',Vl,'faceColor','r','faceAlpha',0.8,'edgeColor','k','lineWidth',0.1);%'faceColor','b');



%close all;
%----------------------- Join node sets
V = [ Vb; Vpp_dx; Vpp_sx;Vpe_dx; Vpe_sx; Vl ];
%F2=[F2+size(V1,1)];
Epp_dx = [ Epp_dx+size(Vb,1) ] ;
Epp_sx = [ Epp_sx+size(Vb,1)+size(Vpp_dx,1) ] ;
Epe_dx = [ Epe_dx+size(Vb,1)+size(Vpp_dx,1)+size(Vpp_sx,1)] ;
Epe_sx = [ Epe_sx+size(Vb,1)+size(Vpp_dx,1)+size(Vpp_sx,1)+size(Vpe_dx,1) ] ;
El = [ El+size(Vb,1)+size(Vpp_dx,1)+size(Vpp_sx,1)+size(Vpe_dx,1)+size(Vpe_sx,1) ] ; 

E = [ Eb; Epp_dx; Epp_sx;Epe_dx; Epe_sx; El];
%----------- total alphashape
shp = alphaShape(V(:,1),V(:,2),V(:,3),4);
%plot(shp)
[facets,nodes] = boundaryFacets(shp);
%figure
%patch('Faces',facets,'Vertices',nodes,'FaceColor','b','lineWidth',0.2,'lineWidth',0.2,'edgeColor','k');

% ----------- improving the surfaces
%find 3points connected
[N] = numConnect(facets,nodes);
logicThree=N==3;
faceBoundMarker1=2;  

[Ft,Vt,~,L]=triSurfRemoveThreeConnect(facets,nodes,[]);
C=double(L);

%surface smoothening
%plot3(Vt(:,1),Vt(:,2),Vt(:,3),'.b');
cPar.n=25;      %check the best value
cPar.Method='HC';
[Vt]=patchSmooth(Ft,Vt,[],cPar);
%Ft = facets;
%plot3(Vt(:,1),Vt(:,2),Vt(:,3),'.b'),hold on;
%patch('Faces',Ft,'Vertices',Vt,'lineWidth',edgeWidth,'edgeColor','c');


% -------- Internal volume for the cancellous bone
% scale the coordinates and get the inner vertebral body, in order to have
% a cortical shell and a cancellous core
V_inn = Vb*0.9;
shp = alphaShape(V_inn(:,1),V_inn(:,2),V_inn(:,3),4);
%plot(shp)
[facets_inn,nodes_inn] = boundaryFacets(shp);
faceBoundMarker2=3;  


%remove3connected points
[F_inn,V_inn,~,L] = triSurfRemoveThreeConnect(facets_inn,nodes_inn,[]);
C=double(L);

%surface smoothening
%plot3(Vt(:,1),Vt(:,2),Vt(:,3),'.b');
cPar.n=25;      %check the best value
cPar.Method='HC';
[V_inn]=patchSmooth(F_inn,V_inn,[],cPar);
%Ft = facets;
%plot3(Vt(:,1),Vt(:,2),Vt(:,3),'.b'),hold on;
%patch('Faces',F_inn,'Vertices',V_inn,'lineWidth',edgeWidth,'edgeColor','c');


V_mix=[Vt;V_inn]; %Joining nodes
F_mix=[Ft;F_inn+size(Vt,1)]; %Joining faces
faceBoundaryMarker=[faceBoundMarker1*ones(size(Ft,1),1); faceBoundMarker2*ones(size(F_inn,1),1)]; %Create boundary markers for faces
if(j==1) %plot only the first
    %plotting surface models
    hf=cFigure;
    title('Surface model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    patch('Faces',Ft,'Vertices',Vt,'FaceColor','r','FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    colormap(autumn(2));
    camlight headlight;
    set(gca,'FontSize',fontSize);
   axis tight;  axis equal;  %grid on;
end

outputpathstlName=fullfile(defaultFolder,'data','output','stl'); 

%write the stl
if stl==1
    outputstlName=fullfile(outputpathstlName,sprintf('vert%d.stl',j));
    %name = sprintf('vert%d',j)
    stlStruct.solidNames={outputstlName};
    stlStruct.solidVertices={Vt};
    stlStruct.solidFaces={Ft};
    stlStruct.solidNormals={[]};
    fileName=fullfile(outputstlName);%'./output/stl/vert%d.stl',j));

    export_STL_txt(fileName,stlStruct);
end

%========================== tetgen  ====================================
fprintf('------------ writing the tetgen file---------------\n');
%filePath=mfilename('fullpath');
%savePath=fullfile(fileparts(filePath));%,'data','temp');
modelName = sprintf('vertebrae_%d',j);
%modelName=('vertebrae');

faceBoundaryMarker=[faceBoundMarker1*ones(size(Ft,1),1); faceBoundMarker2*ones(size(F_inn,1),1)]; %Create boundary markers for faces

V_regions=[mean(Vl(:,1)) mean(Vl(:,2)) mean(Vl(:,3));mean(Vb(:,1)) mean(Vb(:,2)) mean(Vb(:,3)) ]; % to define for both the regions
%plot_matrix(V_regions,'Or');

y1=min(L_body(:,2)) ;
y2 =max(L_lam(:,2)) ;
y_h = y1+(y2-y1)/2 ;
%x1=min(L_lam(:,1));x2 =max(L_lam(:,1));x_h = x1+(x2-x1)/2
z=mean(L_body(:,3));
V_holes=[0,y_h,z];                                   %define the hole

%For each region the mesh density parameter can be specified
[v1]=tetVolMeanEst(Ft,Vt); %Estimate volume of ideal tetrahedron
%regionA=[v1]; % Regional mesh parameters

%For each region the mesh density parameter can be specified
[v2]=tetVolMeanEst(F_mix,V_mix); %Estimate volume of ideal tetrahedron
regionA=[v1*2 v2*4];



stringOpt='-pq2.0AaYQ';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F_mix;
inputStruct.Nodes=V_mix;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%[meshOutput]=runTetGen(inputStruct); %Run tetGen
[meshOutput]=runTetGen(inputStruct);

FT=meshOutput.faces;
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
VT=meshOutput.nodes;
C=meshOutput.faceMaterialID;
E=meshOutput.elements;
elementMaterialIndices=meshOutput.elementMaterialID;

Lmes(j).FT = FT;
Lmes(j).Fb = Fb; 
Lmes(j).Cb = Cb; 
Lmes(j).VT = VT;
Lmes(j).C = C;
Lmes(j).E = E; 
Lmes(j).elementMaterialIndices = elementMaterialIndices;
%plot3(V(:,1),V(:,2),V(:,3),'.r'),hold on;
%figure
%patch('Faces',FT,'Vertices',VT,'FaceColor','b','faceAlpha',1,'edgeColor','k','lineWidth',0.1);%'flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);

%  plot and internal side
%{
hf1=figuremax(figColor,figColorDef);
subplot(1,2,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',FT,'Vertices',VT,'FaceColor','flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);


%hf3=figuremax(figColor,figColorDef);

%subplot(1,3,3);
%Selecting half of the model to see interior
%{
Y=VT(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),C(L));
%{
title('Cut view of solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',VT,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor);
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;
%}
%}

%}
end

%name = sprintf('./output/mat/Lmes_v2.mat');

fid = fopen(outputName,'w');
save(outputName, 'Lmes');
fclose(fid);
%---------------end vertebrae pre-processing


%% 
%
% 
% 
% _*LMG*_ 
% <https://celavecchia.github.io/LMG/>
% 
% _Carolina Eleonora Lavecchia_, <lavecchia.carolina@gmail.com>

