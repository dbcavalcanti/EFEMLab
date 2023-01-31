%% EFEMDraw class
%
% This class implements methods to plot graphical results
% from the EFEM analysis
%
%% Author
% Danilo Cavalcanti
%
%% Class definition
classdef EFEMDraw < handle
    %% Public properties
    properties
        model = []; % handle to an object of the EFEMmodel class
    end
    
    %% Class (constant) properties
    properties (Constant)
        supsize_fac = 0.035; % factor for support size
        loadsize_fac = 0.3; % factor for maximum load size in relation to
                             % half vertical window size
        minloadsize_fac = 0.012 % factor for minimum load size
        arrowsize_fac = 0.025; % factor for load arrow size
        loadstep_fac = 0.05; % factor for load step
        dimlineshift_fac = 0.85; % factor for down shift of dimension lines
                                 % with respect to half Y size
        dimlinetick_fac = 0.15; % factor for dimension line tick sizeEFEM
                                % with respect to half Y size
        maxdisplsize_fac = 0.60; % factor for maximum transversal size in
                                 % relation to half vertical window size
        inflectpt_fac = 0.005; % factor for size of inflection point disk
        rotmeter_fac = 0.85; % factor for rotation meter size in
                             % relation to half vertical window size
        picktol_fac = 0.01; % factor for picking a point
        minmemblen_fac = 0.05; % factor for minimum member length
        ValidSupInsertion = 1; % status for valid support insertion
        BeamLineNotFound = 2; % status for beam line not found for support insertion
        SupInsertionNotValid = 3; % status for not valid position for support insertion
        SupDelMinNumSup = 1; % status for minimum number of internal supports
        ValidSupDeletion = 2; % status for valid support deletion
        SupDelNotFound = 3; % status for support not found for deletion
        MembLoadFound = 1; % status for pick member load found
        MembLoadNotFound = 2 % status for pick member not found
        SupMoveFound = 1; % status for support found for moving
        SupMoveNotFound = 2; % status for support not found for moving
    end
    
    %% Private properties
    properties (Access = private)
        pickmember = 0; % current member load picked
        picksup = 0; % current support picked
        orig_suppos = 0; % original moving support position
        hnd_draft = []; % handle to draft graphics object being displayed
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function draw = EFEMDraw(model)
            if (nargin > 0)
                draw.model = model;
            end
        end
    end
    
    %% Class (static) auxiliary functions
    methods (Static)
        %------------------------------------------------------------------
        % Plots a square with defined center coordinates, side length and
        % color.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  l: side length
        %  c: color (RGB vector)
        function square(cnv,x,y,l,c)
            X = [x - l/2, x + l/2, x + l/2, x - l/2];
            Y = [y - l/2, y - l/2, y + l/2, y + l/2];
            fill(cnv, X, Y, c);
        end
        
        %------------------------------------------------------------------
        % Plots a draft version of a square with defined center coordinates,
        % side length and color.
        % It draws a hollow square and sets its parent property as the
        % given handle to the group of graphics objects.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  hnd: handle to group of graphics objects
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  l: side length
        %  c: color (RGB vector)
        function draftSquare(cnv,hnd,x,y,l,c)
            X = [x - l/2, x + l/2, x + l/2, x - l/2, x - l/2];
            Y = [y - l/2, y - l/2, y + l/2, y + l/2, y - l/2];
            plot(cnv, X, Y, 'color', c, 'Parent', hnd);
        end
        
        %------------------------------------------------------------------
        % Plots a triangle with defined top coordinates, height, base,
        % orientation, and color.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: triangle height
        %  b: triangle base
        %  ang: angle (in radian) between the axis of symmetry and the
        %       horizontal direction (counterclockwise) - 0 rad when
        %       triangle is pointing left
        %  c: color (RGB vector)
        function triangle(x,y,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + h * cx + b/2 * cy, x + h * cx - b/2 * cy];
            Y = [y, y + h * cy - b/2 * cx, y + h * cy + b/2 * cx];
            fill(X, Y, c);
        end
        
        %------------------------------------------------------------------
        % Plots a draft version of a triangle with defined top coordinates,
        % height, base, orientation, and color.
        % It draws a hollow triangle and sets its parent property as the
        % given handle to the group of graphics objects.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  hnd: handle to group of graphics objects
        %  x: top coordinate on the X axis
        %  y: top coordinate on the Y axis
        %  h: triangle height
        %  b: triangle base
        %  ang: angle (in radian) between the axis of symmetry and the
        %       horizontal direction (counterclockwise) - 0 rad when
        %       triangle is pointing left
        %  c: color (RGB vector)
        function draftTriangle(hnd,x,y,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + h * cx + b/2 * cy, x + h * cx - b/2 * cy, x];
            Y = [y, y + h * cy - b/2 * cx, y + h * cy + b/2 * cx, y];
            plot(X, Y, 'color', c, 'Parent', hnd);
        end
        
        %------------------------------------------------------------------
        % Plots a circle with defined center coordinates, radius and color.
        % This method is used to draw hinges on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: circle radius
        %  c: color (RGB vector)
        function circle(cnv,x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            plot(cnv, xcirc, ycirc, 'color', c);
        end
        
        %------------------------------------------------------------------
        % Plots a circle disk with defined center coordinates, radius and
        % color. The circle is filled with the given color
        % This method is used to inflection on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: center coordinate on the X axis
        %  y: center coordinate on the Y axis
        %  r: circle radius
        %  c: color (RGB vector)
        function disk(cnv,x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            fill(cnv, xcirc, ycirc, c);
        end
        
        %------------------------------------------------------------------
        % Plots an arrow with defined beggining coordinates, length,
        % arrowhead height, arrowhead base, orientation, and color.
        % This method is used to draw load symbols on 2D models.
        % Input arguments:
        %  cnv: graphics context (axes)
        %  x: beggining coordinate on the X axis
        %  y: beggining coordinate on the Y axis
        %  l: arrow length
        %  h: arrowhead height
        %  b: arrowhead base
        %  ang: pointing direction (angle in radian with the horizontal
        %       direction - counterclockwise) - 0 rad when pointing left
        %  c: color (RGB vector)
        function arrow2D(cnv,x,y,l,h,b,ang,c)
            cx = cos(ang);
            cy = sin(ang);

            X = [x, x + l * cx];
            Y = [y, y + l * cy];
            line(cnv, X, Y, 'Color', c,'LineWidth',2);
            EFEMDraw.triangle(x, y, h, b, ang, c);
        end

        %------------------------------------------------------------------
        % Snap a value to the closest step value.
        function snap_val = snapToStepValue(val,step)
            fp = val / step;   % "fraction" part
            ip = floor(fp);    % integer part
            fp = fp - ip;
            if fp > 0.5
                snap_val = (ip + 1.0) * step;
            elseif fp < -0.5
                snap_val = (ip - 1.0) * step;
            else
                snap_val = ip * step;
            end
        end
    end

    %% Protect methods
    methods (Access = public) % Access from methods in subclasses
        
        function max_load = getMaxLoad(draw)
            max_load = max(draw.model.F);
        end

        function elements(draw)
            for el = 1:draw.model.nelem
                res = draw.model.element(el).type.result;
                patch('Faces',res.faces,...
                    'Vertices',res.vertices,...
                    'FaceVertexCData',res.vertexData,...
                    'FaceColor',res.faceColor,...
                    'LineWidth',res.edgesThickness,...
                    'Marker',res.markerType,...
                    'MarkerFaceColor',res.markerFaceColor); 
                colormap(res.colormapType);
                if isa(draw.model.element(el).type,'EnrichedElement')
                    res = draw.model.element(el).type.fracture.result;
                    patch('Faces',res.faces,...
                        'Vertices',res.vertices,...
                        'FaceColor','white',...
                        'EdgeColor','k',...
                        'EdgeAlpha',0.5,...
                        'LineStyle','--',...
                        'LineWidth',res.edgesThickness,...
                        'Marker',res.markerType,...
                        'MarkerFaceColor',res.markerFaceColor); 
                end

            end
        end

        function fractures(draw)
            for el = 1:size(draw.model.FRACT,1)
                x = [draw.model.NODE_D(draw.model.FRACT(el,:),1)];
                y = [draw.model.NODE_D(draw.model.FRACT(el,:),2)];
                p = plot(x,y,'k--o','MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1);
                p.Color(4) = 0.3;
            end
        end
    end
    
    %% Public methods
    methods        
        %------------------------------------------------------------------
        % Returns the bounding box (x and y limits) of a continuous beam
        % model. The returned box has xmin = 0, xmax = totalLen,
        % ymin = -totalLen * 0.05, ymax = +totalLen * 0.05, in which
        % totalLen is the length of the entire beam model.
        % The y limits are fictitious. They are equal in module and 
        % equal to a small percentage of the total length to force 
        % the adjustment of the box in the y direction, keeping y = 0
        % in the center of the canvas.
        function bbox = EFEMBoundBox(draw)
            minX = Inf; minY = Inf; maxX = -Inf; maxY = -Inf;
            for el = 1:draw.model.nelem
                res = draw.model.element(el).type.result;
                minX = min(minX,min(res.vertices(:,1)));
                minY = min(minY,min(res.vertices(:,2)));
                maxX = max(maxX,max(res.vertices(:,1)));
                maxY = max(maxY,max(res.vertices(:,2)));
            end
            tolx = 0.10*abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1)));
            toly = 0.10*abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1)));
            bbox(1) = minX - tolx;
            bbox(2) = minY - toly;
            bbox(3) = maxX + tolx;
            bbox(4) = maxY + toly;
        end
        
        %------------------------------------------------------------------
        % Draws a continuous beam model with applied loads.
        % Input:
        % - cnv: graphics context (owning canvas)
        function mesh(draw)

            figure
            hold on, box on, grid on
            axis equal
            
            % Draw continuous mesh
            draw.elements();

            % Draw fractures
            draw.fractures();

            % Draw supports
            dim = max(abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1))),abs(max(draw.model.NODE(:,1)) - min(draw.model.NODE(:,1))));
            supsize = dim * EFEMDraw.supsize_fac;
            for i = 1:draw.model.nnodes
                if (draw.model.SUPP(i,1)==1 || draw.model.SUPP(i,2)==1)
                    posx = draw.model.NODE(i,1);
                    posy = draw.model.NODE(i,2);
                    EFEMDraw.triangle(posx,posy,supsize,supsize,-pi/2,[0 0 1]);
                end
            end

            % Get the bounding box
            bbox = draw.EFEMBoundBox();
            xlim([bbox(1) bbox(3)]);
            ylim([bbox(2) bbox(4)]);
            
        end
        
        %------------------------------------------------------------------
        % Cleans data structure of a EFEMmodel object.
        function draw = clean(draw)
            draw.model = [];
        end
    end
end