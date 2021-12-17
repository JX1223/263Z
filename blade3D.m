function blade3D(a, c, const)
    fp = ['.' filesep 'figures' filesep];
    % Colors
    color_skin = [59, 126, 161]./255; % Founders rock
    %color_skin = [0 50/255 98/255]; % Berkeley blue
    color_section = [253/255, 181/255, 21/255]; % California gold

    % Open file
    fID = fopen("airfoil.scad","w+");
    
    % Calculate radius
    dr = (const.r_t - const.r_r)/(numel(a)); % Width of each blade element
    r = linspace(const.r_r, const.r_t, numel(a)+1)'; 
    r = r(1:end-1) + dr/2; % Distance to 'midpoint' of each blade element
    % Calculate angle phi
    V_t = const.omega * r;
    phi = atan(const.V_h./V_t);

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
    f1 = figure;
    view(3);
    hold on;
    xlabel("x");
    ylabel("r [m]");
    zlabel("y");
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
        if i == 1 || mod(i, 5) == 0
            pg = plot(section, "Parent", t(i));
            pg.FaceColor = color_section;
            pg.FaceAlpha = 0;
            pg.EdgeColor = color_section;
            pg.LineWidth = 1.5;
        end
    end
    fprintf(fID,"];");
    fclose(fID);
    % Generate stl file
    
    if ismac
        status = system("/Applications/OpenSCAD.app/Contents/MacOS/OpenSCAD -q -o blade.stl generate_airfoil.scad");
    elseif ispc
        cmd = sprintf("cd C:\\Program Files\\OpenSCAD\\ & openscad -q -o %s\\blade.stl %s\\generate_airfoil.scad", pwd, pwd);
        status = system(cmd);
    end
    % Import stl mesh
    blade_TR = stlread('blade.stl');
    p = trimesh(blade_TR,'FaceColor',color_skin,'EdgeColor','none','FaceAlpha',0.4);
    rotate(p,[1 0 0], -90, [0, 0, 0]);
    rotate(p,[0 1 0], 180, [0, 0, 0]);
    v = get(p, 'vertices');
    %v(:,1) = v(:,1) + 0;
    %set(p,'vertices',v);
    %ylim([0, 10]);
    axis equal;
    xlim([min(v(:,1)), max(v(:,1))]);
    ylim([1, 10]);
    grid on;
    hold off;
    %view([45 - 180, 45]);
    view([-155, 16]);
    exportgraphics(f1, strcat(fp, '3D.eps'));
end