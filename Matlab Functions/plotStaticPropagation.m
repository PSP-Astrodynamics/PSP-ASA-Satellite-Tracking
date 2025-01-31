function [] = plotStaticPropagation(y)

R = 6378.1363; % km

% 3D plot
plot3(y(:,1),y(:,2),y(:,3), 'r')
hold on
%plot3(x_moon(:,1),x_moon(:,2),x_moon(:,3), 'w')
%plot3(x_sun(:,1),x_sun(:,2),x_sun(:,3), 'y')
%plot3(x_sat(:,1),x_sat(:,2),x_sat(:,3), 'g')

xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

% visual earth model
% [s1, s2, s3] = sphere; % coordinates of a sphere in cartesian
% surf(R*s1,R*s2,R*s3, 'Edgecolor', 'b'); % plot sphere with radius of planet
earthy(R)
axis equal;
grid on;

set(gca, 'Color','k', 'XColor','w', 'YColor','w', 'ZColor','w')
%set(gca,'Color',[.7 .7 .7])
set(gcf,'color','k','units','normalized','outerposition',[0.5 0 0.5 1]);

% initial location
plot3(y(1,1),y(1,2),y(1,3), 'b*')
axis equal

end