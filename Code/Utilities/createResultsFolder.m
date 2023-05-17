% *************************************************************************
% MAIN AUTHOR: Francisco J. Martin-Vega Francisco J. Martin Vega
% *************************************************************************
% GROUP: Lab 1.3.5., Communications and Signal Processing Lab (ComSP), 
% Telecommunication Research Institute (TELMA), ETSIT, University of Malaga
% *************************************************************************
% DESCRIPTION:
% This script evaluates the performance of power setting
% *************************************************************************

function p = createResultsFolder(p)

[~, message, ~] = mkdir('Results');

strDate = date;
p.folderName = ['Results\' p.resultFolder '_' strDate];
[~, message, ~] = mkdir(p.folderName);
i = 1;
folderNameAux = p.folderName;
while strcmp(message, 'Directory already exists.')
    folderNameAux = [p.folderName ' (' num2str(i) ')'];
    i = i + 1;
    [~, message, ~] = mkdir(folderNameAux);
end
p.folderName = folderNameAux;

[ST, ~] = dbstack; 
copyfile(ST(3).file,[p.folderName '\' 'params_' p.resultFolder '.m']); 

