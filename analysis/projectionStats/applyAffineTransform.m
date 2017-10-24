function transformed = applyAffineTransform(orig,s,target,doPlot)
% Apply output of absor to orig in the hope that it will rotate to look like target
% 
% function transformed =  applyAffineTransform(orig,s,target)    
%
% Inputs
% - orig and target should be 3 by N matrices (as requested by absor)
% - s is the output of absor
% - "target" is optional. It is the target of the rotation os the difference between 
% - the transformed points produced here and target is are residuals.
% - doPlot is true by default
% 
% http://ch.mathworks.com/matlabcentral/fileexchange/26186-absolute-orientation-horn-s-method
%
% More details:
% https://www.mathworks.com/matlabcentral/answers/93554-how-can-i-rotate-a-set-of-points-in-a-plane-by-a-certain-angle-about-an-arbitrary-point
% https://people.cs.clemson.edu/~dhouse/courses/401/notes/affines-matrices.pdf
% 
% Rob Campbell - Basel 2017

if nargin<4 || isempty(doPlot)
    doPlot=true;
end

transformed = zeros(size(orig,1)+1, size(orig,2));

for ii=1:size(orig,2)
    transformed(:,ii) = s.M * [orig(:,ii);1];
end

transformed(end,:)=[]; %Remove w

% Plot
if ~doPlot, return, end

%Cap the number of points to 1000
maxPoints = 1E3;
if size(transformed,2)>maxPoints
    ind = round(linspace(1,size(transformed,2),maxPoints));
else
    ind = 1:size(transformed,2);
end

cla

hold on
for k=1:length(ind)
    ii = ind(k);
    plot3([orig(1,ii), transformed(1,ii)], ...
        [orig(2,ii),transformed(2,ii)], ...
        [orig(3,ii),transformed(3,ii)], ...
        '-', 'Color', [0.73,0.73,1]);
end



plot3(transformed(1,ind), transformed(2,ind), transformed(3,ind),'.r','MarkerSize',15)
plot3(orig(1,ind), orig(2,ind), orig(3,ind),'.k', 'MarkerSize',10)

%overlay the target if present
if nargin>2 && ~isempty(target)
    plot3(target(1,ind), target(2,ind), target(3,ind), 'ob', 'MarkerSize',10)
end


view(3)

hold off

xlabel('x'), ylabel('y'), zlabel('z')
axis equal
grid on
box on
