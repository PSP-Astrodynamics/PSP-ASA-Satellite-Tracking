clc
clear
close all

%% pos
% Load the CSV file
data = readtable('Ephemeris_Position_catsat.csv');

% Extract columns
time = data{:, 1}; % Time column
nx = data{:, 2};    % X position column (in meters)
ny = data{:, 3};    % Y position column (in meters)
nz = data{:, 4};    % Z position column (in meters)

% Skip rows with missing values
validRows = ~isnan(nx) & ~isnan(ny) & ~isnan(nz);

% Convert valid positions to kilometers
nx = nx(validRows) / 1000;
ny = ny(validRows) / 1000;
nz = nz(validRows) / 1000;

r_vec_xyz = [nx,ny,nz];

% Plot the 3D trajectory
% figure;
% plot3(x, y, z, '-o');
% grid on;
% xlabel('X Position (km)');
% ylabel('Y Position (km)');
% zlabel('Z Position (km)');
% title('3D Satellite Position Trajectory');
% hold on
% earthy(6378.1363, 'Earth', 1,[0,0,0])
% axis equal

%% vel
% Load the CSV file
data = readtable('Ephemeris_Velocity_catsat.csv', 'TextType', 'string');

% Extract the velocity columns (X, Y, Z) as strings
vx_str = data{:, 'X'};
vy_str = data{:, 'Y'};
vz_str = data{:, 'Z'};

% Remove 'm/s' from the strings and convert to numbers
vx = str2double(erase(vx_str, 'm/s'));
vy = str2double(erase(vy_str, 'm/s'));
vz = str2double(erase(vz_str, 'm/s'));

% Remove rows with NaN values (corresponding to empty or invalid entries)
validRows = ~isnan(vx) & ~isnan(vy) & ~isnan(vz);
vx = vx(validRows)/1000;
vy = vy(validRows)/1000;
vz = vz(validRows)/1000;

v_vec_xyz = [vx,vy,vz];

%% analys

for n=1:length(r_vec_xyz)
    
    h_vec_xyz(n,:) = cross(r_vec_xyz(n,:),v_vec_xyz(n,:));
    z_vec_xyz(n,:) = [0,0,1];

    nodal_vec_xyz(n,:) = cross(z_vec_xyz(n,:),h_vec_xyz(n,:));

    RAAN(n) = 360-acosd(dot(nodal_vec_xyz(n,:),[1,0,0])/norm(nodal_vec_xyz(n,:)));
end

% Plot the 3D nodal vec
% Extract the components
nx = nodal_vec_xyz(:, 1); % X components
ny = nodal_vec_xyz(:, 2); % Y components
nz = nodal_vec_xyz(:, 3); % Z components

% Create the 3D plot
figure;
hold on;
grid on;

% Plot a vector from the origin to each point
for i = 1:size(nodal_vec_xyz, 1)
    plot3([0, nx(i)], [0, ny(i)], [0, nz(i)], '-o'); % Plot line and endpoint
end

% Add labels and title
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
title('3D Nodal Vector for 3-Day Period');
view(3); % Set 3D view
earthy(6378.1363, 'Earth', 1,[0,0,0])
axis equal
hold off;

nodal1_vec_xyz = nodal_vec_xyz(1,:);
n1 = norm(nodal1_vec_xyz);
nodalend_vec_xyz = nodal_vec_xyz(end,:);
nend = norm(nodalend_vec_xyz);

DRAAN = RAAN(end)-RAAN(1); % deg

t = juliandate(time(validRows));
p = polyfit(t,RAAN,1);
RAANfit = polyval(p,t);

figure
plot(time(validRows),RAAN,'bo')
grid on
title('Right Ascension of the Ascending Node vs Time')
xlabel('time')
ylabel('RAAN (\Omega) [deg]')
hold on
plot(time(validRows),RAANfit, '--')
legend('RAAN measured','Linear fit')
