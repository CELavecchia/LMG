function install_LMG

%% Add GIBBON paths
disp('---INSTALLING THE LMG TOOLBOX ---');
disp('-> Adding LMG paths');

toolboxPath=fileparts(mfilename('fullpath')); %Get the LMG path
addpath(fullfile(toolboxPath,'lib')); %Add toolbox lib path so toolbox functions used here are known

%%
[pathNames]=getSubPaths(toolboxPath);

%%
ignoreString='.'; %Target for ignoring .git and other hidden folders
for q=1:1:numel(pathNames)
    pathNameNow=pathNames{q};
    if ~any(strfind(pathNameNow,ignoreString))        
        addpath(pathNameNow); %Add path
    end    
end

%% Add 3rd party paths
% disp('-> Adding 3rd party paths');
% 
% prompt = {'FEBio full path to program (e.g. .../bin/FEBio2.lnx64 or ...\bin\FEBio2.exe):','export_fig path:'};
% dlg_title = 'Path definitions (leave empty if not used)';
% 
% FEBioPath=getFEBioPath;
% 
% exportFigPath=fileparts(which('export_fig'));
% defaultOptions = {FEBioPath,exportFigPath};
% 
% s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
% 
% Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
% 
% if ~isempty(Q)
%     if ~isempty(Q{1})
%         setFEBioPath(Q{1}); %Set FEBio path in config file                
%     end        
%     if ~isempty(Q{2})
%         addpath(Q{2}); %Add export_fig to the path             
%     end        
% end

%% Saving path definitions
disp('-> Saving path definitions');
savepath;

%% Integrating help/documentations

disp('-> Integrating help');
createHelpDemoDocumentation;

disp('------- FINISHED -------');
