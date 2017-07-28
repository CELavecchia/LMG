%clear; clc; close all
function meshStruct = Disc_mesh_KMM_gui(EP); 

%% plot settings
fontSize=15;
faceAlpha1=0.2;
markerSize=25;
lineWidth=2;
plotColors=gjet(10); %Plot colors and colormap

%% Control parameters
n_points = 4*16; %number of points on the perimeter, should be multiple of 4!
numHexElements_zDir = 8;
centralShapeOpt=2; %1=square, 2= rectangle
numStepsIntermediateMesh=[]; %Leave empty to base of mesh spacing

interpMethod='pchip'; %Interp method for resampling curves
numDigitsKeepUnique=5; %Used for rounding to merge nodes
radiusFactorInnerMesh=0.5; %Factor to set radius for inner rectangle/square
controlParSmooth.Method='LAP'; %Smoothening method HC or LAP, LAP is more agressive
controlParSmooth.n=25; %Number of surface mesh smoothening steps

%% input

load('./data/input/layers.mat'); 
layers_disc =layers;
size(layers_disc);



layers = layers_disc(:,1:27); %to define according to the VP
nucleus_lay = layers_disc(:,28:end);

vert_step =  8; % elements in z

%Convert to nx3xm 3D array such that 3rd dimension is for the layers
layersRaw=reshape(layers,[size(layers,1),3,size(layers,2)/3]);

%pre-processing curves
interpFunction = fitting_surfaces_gui(EP);


%layersRaw
%% Evenlysample points for all curves

%{
cFigure; 
title('Resampled input curves','FontSize',fontSize);
hold on; 
%}
%axisGeom;

layersResampled=zeros(n_points,3,size(layersRaw,3));
for q=1:1:size(layersRaw,3)
    V_layer_now=layersRaw(:,:,q);
    
    %plotV(V_layer_now,'k-'); 

    [V_layer_now]=evenlySampleCurve(V_layer_now,n_points,'pchip',1);

    %hp=plotV(V_layer_now,'k.-','MarkerSize',markerSize); 
    %set(hp,'Color',plotColors(q,:));
    
    layersResampled(:,:,q)=V_layer_now;
    drawnow; 
end

%% Build surface from layers
% (largely copied over from polyLoftLinear)

%Coordinate matrices
X=squeeze(layersResampled(:,1,:))';
Y=squeeze(layersResampled(:,2,:))';
Z=squeeze(layersResampled(:,3,:))';

c=(1:1:size(Z,1))';
C_layers=c(:,ones(1,size(Z,2)));

%Create quad patch data
[F_layers,V_layers,C_layers] = surf2patch(X,Y,Z,C_layers);

I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
F_sub=sub2ind(size(Z),I,J);
F_layers=[F_layers;F_sub];
[C_layers]=vertexToFaceMeasure(F_layers,C_layers);
C_layers(end-size(F_sub,1):end,:)=C_layers(end-size(F_sub,1):end,:)+0.5;

C_layers=round(C_layers);

V_inner=V_layers(size(layersRaw,3):size(layersRaw,3):end,:); %The inner curve

%% Visualising
%{
cFigure; hold on;
title('Meshed top without anulus','FontSize',fontSize);

gpatch(F_layers,V_layers,C_layers,'k',0.5);
plotV(V_inner,'b.-','MarkerSize',markerSize,'lineWidth',lineWidth)

%axisGeom(gca,fontSize);
colormap(plotColors);
camlight headlight;
drawnow;
%}
%% Build interpolation function for z-coordinates of annulus
interpFunc = scatteredInterpolant(V_layers(:,[1 2]),V_layers(:,3),'natural');

%% Building anulus mesh

[thetaInner,r_inner] = cart2pol(V_inner(:,1),V_inner(:,2));

% Creating central regular quad mesh
nElements=(n_points/4)+1;
[X_centralMesh,Y_centralMesh]=meshgrid(linspace(-1,1,nElements));
[F_centralMesh,V_centralMesh] = surf2patch(X_centralMesh,Y_centralMesh,zeros(size(X_centralMesh)));
V_centralMesh=V_centralMesh(:,1:2);

% Scaling radius
[thetaMesh,radiusMesh]=cart2pol(V_centralMesh(:,1),V_centralMesh(:,2));

radiusMesh=radiusFactorInnerMesh*(1/2)*sqrt(2)*radiusMesh;

switch centralShapeOpt
    case 1 %square
        radiusMesh=mean(r_inner)*radiusMesh;
        [V_centralMesh(:,1),V_centralMesh(:,2)]=pol2cart(thetaMesh,radiusMesh);
    case 2 %rectangle ellipse fit based
        [V_centralMesh(:,1),V_centralMesh(:,2)]=pol2cart(thetaMesh,radiusMesh);
        V_centralMesh(:,3)=0;
        
        [A] = ellipseFit(V_inner(:,[1 2])); %Fit ellipse to inner curve
        x0=A(1);
        y0=A(2);
        V_centralMesh(:,1)=A(3).*V_centralMesh(:,1);
        V_centralMesh(:,2)=A(4).*V_centralMesh(:,2);
        [R,~]=euler2DCM([0 0 -A(5)]);
        V_centralMesh=(R*V_centralMesh')';
        V_centralMesh(:,1)=V_centralMesh(:,1)+x0;
        V_centralMesh(:,2)=V_centralMesh(:,2)+y0;
end
V_centralMesh(:,3)=interpFunc(V_centralMesh(:,[1 2])); 

%gpatch(F_centralMesh,V_centralMesh,'r','k',0.5);

Eb=patchBoundary(F_centralMesh,V_centralMesh);
indOuter=edgeListToCurve(Eb);
indOuter=indOuter(1:end-1);

[~,indMin]=min(abs(thetaMesh(indOuter)-thetaInner(1)));
if indMin>1
    indOuter=[indOuter(indMin:end) indOuter(1:indMin-1)];
end
%plotV(V_centralMesh(indOuter,:),'b.-','MarkerSize',markerSize,'lineWidth',lineWidth)
% plotV(V_centralMesh(indOuter(indMin),:),'r.','MarkerSize',markerSize,'lineWidth',lineWidth)

%%

controlParLoftLinear.numSteps=numStepsIntermediateMesh;
controlParLoftLinear.closeLoopOpt=1;
controlParLoftLinear.patchType='quad';
[F_interMediate,V_interMediate]=polyLoftLinear(V_inner,V_centralMesh(indOuter,:),controlParLoftLinear);

V_interMediate(:,3)=interpFunc(V_interMediate(:,[1 2])); 

%gpatch(F_interMediate,V_interMediate,'r','k',0.5);

%% Join and merge surface data sets

%Color labels for face sets, change as desired
%C_layers=ones(size(F_layers,1),1);
C_interMediate=10*ones(size(F_interMediate,1),1);
C_centralMesh=11*ones(size(F_centralMesh,1),1);

%V_intermediate and V_centralMesh are parts of the nucleus
[Fq,Vq,Cq]=joinElementSets({F_layers,F_interMediate,F_centralMesh},{V_layers,V_interMediate,V_centralMesh},{C_layers,C_interMediate,C_centralMesh});

%merge node sets
[~,indSelect,indFix]=unique(pround(Vq,numDigitsKeepUnique),'rows'); 
Vq=Vq(indSelect,:);
Fq=indFix(Fq);
Eb=patchBoundary(Fq,Vq);

%%
%{
cFigure; hold on;
title('Merged surface geometry','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

gpatch(Fq,Vq,Cq,'k',1);
plotV(Vq(Eb(:),:),'r.','MarkerSize',markerSize,'lineWidth',lineWidth)

%axisGeom(gca,fontSize);
colormap(plotColors);
camlight headlight;
drawnow;
%}
%% Smoothening mesh and reinterpolating (not required but improves mesh)

if controlParSmooth.n>0
    logicFacesHoldOn=(Cq==1);
    indHoldOn=Fq(logicFacesHoldOn,:);
    indHoldOn=unique(indHoldOn(:));
    
    logicFacesSmooth=~logicFacesHoldOn;
    indSmooth=Fq(logicFacesSmooth,:);
    indSmooth=unique(indSmooth(:));
       
    controlParSmooth.RigidConstraints=indHoldOn;    
    [Vq]=patchSmooth(Fq,Vq,[],controlParSmooth); %Smooth
    Vq(indSmooth,3)=interpFunc(Vq(indSmooth,[1 2])); %Fix z coordinate again
end
%%
%{
cFigure; hold on;
title('Smoothened annulus','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

gpatch(Fq,Vq,Cq,'k',1);
% patchNormPlot(F,V); %Plot face normals to ensure they are upward
plotV(Vq(Eb(:),:),'r.','MarkerSize',markerSize,'lineWidth',lineWidth)

%axisGeom(gca,fontSize);
colormap(plotColors);
camlight headlight;
drawnow;
%}
%% Creating hexahedral mesh

offSetDir=[0 0 -1]; %use 1 for facenormaldir, -1 for opposite or provide other directions
[E,V]=quadThick(Fq,Vq,[0 0 -1],vert_step,numHexElements_zDir);
C=repmat(Cq,numHexElements_zDir,1); %Material indices (outer, intermediate, grid)

[FE,CF]=element2patch(E,C,'hex8'); %Set of all element faces
[indFree]=tesBoundary(FE,V); %Indices for boundary faces
Fb=FE(indFree,:); %The boundary faces

% Create faceBoundaryMarkers based on normals (prone to error if surface is very curved so do check)
[N]=patchNormal(Fb,V); %N.B. Change of convention changes meaning of front, top etc.
faceBoundaryMarker=3*ones(size(Fb,1),1); %
faceBoundaryMarker(N(:,3)>0.5)=1; %Top
faceBoundaryMarker(N(:,3)<-0.5)=2; %Bottom



%force the upper and inferior surface





%%
%{
cFigure; 
subplot(1,2,1); 
hold on;
title('Boundary of model with faces colored','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

gpatch(Fb,V,faceBoundaryMarker,'k',1);

%axisGeom(gca,fontSize);
colormap(plotColors);
camlight headlight;

subplot(1,2,2); 
hold on;
title('Cut view of hexahedral mesh','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);

%Create cut view
Y=V(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),C(L),'hex8');
gpatch(Fs,V,Cs,'k',1);
gpatch(Fb,V,'b','none',faceAlpha1);
%axisGeom(gca,fontSize);
camlight headlight;
colormap(plotColors);
drawnow;
%}
%% Output

meshStruct.E=E; %Hex elements
meshStruct.C=C; %Element colors
meshStruct.V=V;
meshStruct.FE=FE;
meshStruct.indFree=indFree;
meshStruct.Fb=Fb;
meshStruct.faceBoundaryMarker=faceBoundaryMarker;
%}
%%
