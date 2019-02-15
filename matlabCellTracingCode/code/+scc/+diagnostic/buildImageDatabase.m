function dataStruct=buildImageDatabase(X,suppressImageCreation)
% Make an HTML page that shows all our traced cells
%
% function dataStruct=scc.diagnostic.buildImageDatabase(X,suppressImageCreation)
%
%
% Purpose
% Call from project root dir to make the webpage
%
% Inputs
% X is the xylem data object
% supressImageCreation is false by default

if nargin<2 || isempty(supressImageCreation)
    suppressImageCreation=false;
end

if ~suppressImageCreation
    scc.diagnostic.buildImages(X)
end

scc.diagnostic.buildHTML(X)
