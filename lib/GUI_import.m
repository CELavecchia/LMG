%% GUI for the LMG
% - import sex, age and height for the correlation analysis;
% - import the curvature of the lumbar spine
%clc; clear all; close all;

function [varagout] = GUI_import
%default path
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
outputpathName=fullfile(defaultFolder,'data','input');

dlg_title = 'LMG input: select model type';

% 1. select the type of model
prompt = {'Select the type of model: 0 for average dimensions, 1 for subject specific measurements '};
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if(str2num(answer{1,1})==0)

    dlg_title = 'LMG input: dimensions based on correlation analysis';

    % 1. select the type of model
    prompt = {'Sex of the patient (m/f)','Enter the age [years]',...
        'Enter the height [cm]','Enter the lumbar curvature [degrees]'};
    num_lines = 1;
    defaultans = {'m','30','180','43.49'};
    answer2 = inputdlg(prompt,dlg_title,num_lines,defaultans);
   
    varagout.type = str2num(answer{1,1});
    varagout.sex = answer2{1,1};
    varagout.age = str2num(answer2{2,1});
    varagout.height = str2num(answer2{3,1});
    varagout.alpha = str2num(answer2{4,1});
    
    %varagout = data;
elseif str2num(answer{1,1})==1 

    dlg_title = 'LMG input';

    % 1. select the type of model
    prompt = {'Sex of the patient (m/f)','Enter the age [years]',...
        'Enter the height [cm]','Enter the lumbar curvature [degrees]','complete path for the xlsx file'};
    num_lines = 1;
    defaultans = {'m','30','180','43.49','subject_specific_dimensions.xlsx'};
    answer3 = inputdlg(prompt,dlg_title,num_lines,defaultans);
    fileName= answer3{5,1};
        
    M=xlsread(sprintf('C:/Users/lavecchc/Documents/LMG/data/input/%s',fileName));
    
    dimensions.type = str2num(answer{1,1});
    dimensions.sex = answer3{1,1};
    dimensions.age = answer3{2,1};
    dimensions.height = answer3{3,1};
    dimensions.alpha = answer3{4,1};
    %read information from the xlsx file
    dimensions.EPWu_half = M(:,1)/2 ;%using the half in the script
    dimensions.EPDu= M(:,2) ;
    dimensions.EPWi= M(:,3) ;
    dimensions.EPDi= M(:,4);
    dimensions.hL=M(:,5); 
    dimensions.PDH=M(:,6); 
    dimensions.PDW=M(:,7); 
    dimensions.lam_l=M(:,8); 
    dimensions.PDs=M(:,9); 
    dimensions.PDt=M(:,10); 
    dimensions.sc_d=M(:,11); 
    dimensions.sc_w=M(:,12); 
    dimensions.TP_wu=M(:,13); 
    dimensions.TP_wi=M(:,14); 
    dimensions.alpha=M(1,15); 
    %dimensions.TP_wi=TP_wi;not used

    %IVD
    dimensions.IVD=M(1:4,18);
    dimensions.IVDd=M(1:4,19);
    dimensions.IVDw=M(1:4,20);
    
    varagout = dimensions;
else
     error('Assign the flag for average dimensions (0) or subject-specific dimensions (1)');
end
   

end
