% Parameterization Lumbar vertebra
% inport previously fitted dataset

%clear; clc; close all;
function [ L, mesh_struct_IVD ] = parameterization_gui(VP,dimensions);





%inport the angle' values
%alpha = angles_parameter;                  %alpha(1)=>L1; alpha(5)=>L5
% alpha in radians


load('./input/L.mat');
load('./input/mesh_struct_IVD.mat');





