function [pdem, facets,nodes] = obtain_geometry( L_body)
w = warning ('off','all');
rmpath('./output')
warning(w);

shp = alphaShape(L_body(:,1),L_body(:,2),L_body(:,3),16);
%plot(shp)
[facets,nodes] = boundaryFacets(shp);
[tri,V] = alphaTriangulation(shp);
%tri = triangulateFaces(shp.tri);
%trimesh(tri,V(:,1),V(:,2),V(:,3));


[N]=numConnect(facets,nodes);
logicThree=N==3;
%remove3connected points
[Ft,Vt,~,L]=triSurfRemoveThreeConnect(facets,nodes,[]);
C=double(L);
%surface smoothening
cPar.n=25;      %check the best value
cPar.Method='HC';
[Vt]=patchSmooth(Ft,Vt,[],cPar);

%get the model

model = createpde();

geometryFromMesh(model,Vt',Ft');
msh = generateMesh(model,'GeometricOrder','linear','Hmax',1.5,'Hmin',1.3);
nodes = msh.Nodes;
elems = msh.Elements;
pdem = createpde;
geometryFromMesh(pdem,nodes,elems);

%figure
%generateMesh(model,'GeometricOrder','linear','Hmax',7,'Hmin',1);
%pdegplot(model,'FaceAlpha',1), hold on;


end