function [] = earthy(R)
    npanels = 180;   % Number of globe panels around the equator
    alpha   = 1;     % Globe transparency level, 1 = opaque, through 0 = invisible
    
    % Earth texture image (URL to a suitable texture map)
    image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
    
    % Mean spherical Earth dimensions
    erad    = R;  % Equatorial radius (meters)
    prad    = R;  % Polar radius (meters)
    
    %% Create wireframe globe (without displaying edges)
    [x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
    globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none'); % No edges

    %% Texturemap the globe
    % Load Earth image for texture mapping
    cdata = imread(image_file);
    
    % Apply texture to the globe
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end
