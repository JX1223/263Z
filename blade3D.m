function blade3D(a, c)
    %% Open file
    fID = fopen("airfoil.scad","w+");
    
    %% Constants
    % Air constants
    rho = 1.225; % kg/m^3
    % Problem constants
    omega = 15 * (2*pi/60); % rpm -> rad/s
    r_r = 1; % m, root radius
    r_t = 10; % m, tip radius
    mass = 115;
    g = 9.81;
    T_R = (mass * g)/(4); % N, thrust required per rotor
    A_R = pi*r_t^2; % m^2, area of disk
    V_h = sqrt(T_R/(2*rho*A_R)); % Hover induced velocity
    % Calculate radius
    dr = (r_t - r_r)/(numel(a)); % Width of each blade element
    r = linspace(r_r, r_t, numel(a)+1)'; 
    r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element
    % Calculate angle phi
    V_t = omega * r;
    phi = atan(V_h./V_t);

    % Convert angle of attack alpha to design angle beta
    beta = rad2deg(-1.*(a + phi));

    % Convert chord to scale factors w.r.t chord at root
    s = c./c(1);

    % Load airfoil:
    load("4515_199.mat","coords");
    
    % Build base airfoil, scaled to root chord length
    warning('off', 'MATLAB:polyshape:repairedBySimplify')
    airfoil = polyshape(coords);
    [xc, yc] = centroid(airfoil);
    airfoil = scale(airfoil,c(1),[xc,yc]);
    %offset = airfoil.Vertices;
    %offset = max(abs(offset(:,1)));
    figure;
    view(3);
    hold on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    t = matlab.graphics.primitive.Transform.empty(numel(a), 0);
    fprintf(fID, "airfoil = [");
    for i = 1:numel(a)
        if i == 1
            fprintf(fID,"[");
        else
            fprintf(fID,",[");
        end
        
        % Scale and rotate airfoil section
        section = rotate(airfoil, beta(i), [xc, yc]);
        section = scale(section,s(i),[xc,yc]);
        p = section.Vertices;
        for j = 1:size(p,1)
            if j == 1
                fprintf(fID,"[%0.6f,%0.6f,%0.6f]", p(j,1), p(j,2), r(i));
            else
                fprintf(fID,", [%0.6f,%0.6f,%0.6f]", p(j,1), p(j,2), r(i));
            end
        end
        fprintf(fID,"]");
        r_i = r(i);
        M = makehgtform("yrotate",pi,"xrotate",-pi/2,'translate',[0, 0, r_i]);
        t(i) = hgtransform('Matrix', M);
        pg = plot(section, "Parent", t(i));
        pg.FaceColor = 'w';
        pg.FaceAlpha = 0;
        pg.EdgeColor = 'k';
    end
    fprintf(fID,"];");
    fclose(fID);
    % Generate stl file
    status = system("/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD -q -o blade.stl generate_airfoil.scad");
    % Import stl mesh
    blade_TR = stlread('blade.stl');
    p = trimesh(blade_TR,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    rotate(p,[1 0 0], -90, [0, 0, 0]);
    rotate(p,[0 1 0], 180, [0, 0, 0]);
    %v = get(p, 'vertices');
    %v(:,1) = v(:,1) + 0;
    %set(p,'vertices',v);
    %ylim([0, 10]);
    axis equal;
    xlim([-2, 1]);
    grid on;
    hold off;
    campos([-24.8779   50.8778   15.0077]);

end