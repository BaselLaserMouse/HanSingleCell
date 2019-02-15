function varargout=makeLegend(addV1)
% Makes a legend with boxes and names color-coded according to brainAreaNames.visualAreas
%
% function H = brainAreaNames.makeLegend(addV1)
%
%
% Usage example:
% H = brainAreaNames.visualAreas
%
%
% Rob Campbell - Basel 2017

if nargin==0
    addV1=false;
end


[areas,colorMap,abbreviations] = brainAreaNames.visualAreas(addV1);

boxSize = [0.15,0.05];

xO = 0.05;
yO = 0.06;


clf
hold on
for ii=1:length(abbreviations)-1 %Skip RS
    thisArea = strrep(areas{ii},' ','_');

    x = [xO, xO+boxSize(1), xO+boxSize(1), xO];
    y = [yO*(ii-1), yO*(ii-1), boxSize(2)+yO*(ii-1), boxSize(2)+yO*(ii-1)];


    H(ii).ptch=patch(x, y, colorMap.areas(thisArea) );
    H(ii).txt=writeText(x(1)+xO/2, (y(1)+y(3))/2, abbreviations{ii});

end


axis off

xlim([0,1])
if nargout>0
    varargout{1}=H;
end


function h=writeText(x,y,txt)
    fSize=20;
    h=text(x,y,txt,'fontSize',fSize,'fontWeight','bold','color',[1,1,1]*0.8);
    h=text(x-0.0015,y+0.0015,txt,'fontSize',fSize,'fontWeight','bold','color','k');