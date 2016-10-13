function hh=plot_h2(dat)

load('flatHelmetSK.mat');
%figure;
%set(gcf,'position',[3          54        1493        1066])
hAxes=gca;
hold on;


hold on
hh=patch('Faces',     Faces,'Vertices',  Vertices,'FaceVertexCData', dat,...
    'EdgeColor', 'k','EdgeAlpha', 0.01, 'FaceColor',...
'interp','BackfaceLighting', 'lit','facealpha',1);
hold on
axis('equal')
axis('off')

%colorbar;
set(gca,'fontsize',12,'fontweight','b')
%set(gca,'clim',[-max(abs(dat)) max(abs(dat))])
radii = [Vertices(:,2);Vertices(:,1)];
PlotNoseEars(hAxes, (max(radii)-min(radii))/4, 1);
end