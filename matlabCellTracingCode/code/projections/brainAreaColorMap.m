function C = brainAreaColorMap
    % Return a map structure that stores the assigned color for each brain area.
    %
    %
    % See: colorBrainAreas, returnVisAreas



P=jet(8);
P=[P(1:2:end,:);P(2:2:end,:)]; %Avoid placing too many similar colors near each other.

C=containers.Map;

C('anteromedial_visual_area')=P(1,:);
C('anterolateral_visual_area')=P(2,:);
C('Lateral_visual_area')=P(3,:);
C('rostrolateral_visual_area')=P(4,:);
C('posterolateral_visual_area')=P(5,:);
C('posteromedial_visual_area')=P(6,:);
C('Anterior_area')=P(7,:);
C('Laterointermediate_area')=P(8,:);
