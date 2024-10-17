classdef SakurajimaApp < matlab.apps.AppBase
    % SakurajimaApp: Aplicación para simular la trayectoria de partículas volcánicas
    % Esta app simula y visualiza la trayectoria de partículas expulsadas por el volcán Sakurajima
    % en 2D y 3D, considerando la resistencia del aire.

    properties (Access = public)
        % Componentes de la interfaz de usuario
        UIFigure                 matlab.ui.Figure
        SimulateButton           matlab.ui.control.Button
        MassofprojectileEditField   matlab.ui.control.NumericEditField
        MassofprojectileEditFieldLabel  matlab.ui.control.Label
        DiameterofprojectileEditField  matlab.ui.control.NumericEditField
        DiameterofprojectileEditFieldLabel  matlab.ui.control.Label
        cEditField               matlab.ui.control.NumericEditField
        cEditFieldLabel          matlab.ui.control.Label
        nEditField               matlab.ui.control.NumericEditField
        nEditFieldLabel          matlab.ui.control.Label
        pEditField               matlab.ui.control.NumericEditField
        pEditFieldLabel          matlab.ui.control.Label
        TimestepEditField        matlab.ui.control.NumericEditField
        TimestepEditFieldLabel   matlab.ui.control.Label
        ShotangleEditField       matlab.ui.control.NumericEditField
        ShotangleEditFieldLabel  matlab.ui.control.Label
        InitialvelocityEditField  matlab.ui.control.NumericEditField
        InitialvelocityEditFieldLabel  matlab.ui.control.Label

        UIAxes3D                 matlab.ui.control.UIAxes
        UIAxes2D                 matlab.ui.control.UIAxes
    end

    methods (Access = private)
        function SimulateButtonPushed(app, event)
            % Función ejecutada al presionar el botón Simular
            % Obtiene los valores de entrada y ejecuta las simulaciones 2D y 3D

            % Obtener valores de entrada
            dt = app.TimestepEditField.Value;
            p = app.pEditField.Value;
            n = app.nEditField.Value;
            c = app.cEditField.Value;
            diameter = app.DiameterofprojectileEditField.Value;
            mass = app.MassofprojectileEditField.Value;
            shot_angle = app.ShotangleEditField.Value;
            initial_velocity = app.InitialvelocityEditField.Value;
            
            % Simulación 2D
            [positions2D, positions2D_no_drag, time_values] = run2DSimulation(app, dt, p, n, c, diameter, mass, shot_angle, initial_velocity);
            
            % Simulación 3D
            [positions3D, x, y, z] = run3DSimulation(app, dt, p, n, c, diameter, mass, shot_angle, initial_velocity);
            
            % Graficar resultados
            plot2DResults(app, positions2D, positions2D_no_drag, time_values);
            plot3DResults(app, positions3D, x, y, z);
        end
        
        % Simulación 2D de la trayectoria de la partícula
        % Calcula la trayectoria con y sin resistencia del aire
        function [positions2D, positions2D_no_drag, time_values] = run2DSimulation(app, dt, p, n, c, diameter, mass, shot_angle, initial_velocity)
            g = 9.81;
            total_time = 10; % Tiempo total de simulación (10s)
        
            % Ángulo de disparo a radianes
            shot_angle_rad = deg2rad(shot_angle);
        
            % Componentes iniciales de velocidad (m/s)
            vx = initial_velocity * cos(shot_angle_rad);
            vz = initial_velocity * sin(shot_angle_rad);
        
            % Cálculo de factor de arrastre
            A = pi * (diameter/2)^2; % Área de partícula
            k = 0.5 * p * c * A / mass; % Factor de arrastre (D/m)
        
            % Inicialización de arrays para la simulación
            time = 0:dt:total_time;
            num_steps = length(time);
            positions2D = zeros(2, num_steps);
            positions2D_no_drag = zeros(2, num_steps);
            velocities2D = zeros(2, num_steps);
            
            positions2D(:,1) = [0; 0];
            positions2D_no_drag(:,1) = [0; 0];
            velocities2D(:,1) = [vx; vz];
        
            hit_ground = false;
            hit_ground_no_drag = false;
            final_x_no_drag = 0;
        
            for t = 2:num_steps
                % Con air drag
                if ~hit_ground
                    v = velocities2D(:,t-1);
                    p = positions2D(:,t-1);
                    speed = norm(v);
                    
                    % Calular la aceleración
                    a_drag = -k * speed^(n-1) * v;
                    a = a_drag + [0; -g];
                    
                    % Actualizar velocidad y posición con Método de Euler
                    velocities2D(:,t) = v + a * dt;
                    positions2D(:,t) = p + velocities2D(:,t) * dt;
                    
                    if positions2D(2,t) <= 0
                        % Sacar el tiempo exacto de impacto
                        t_impact = time(t-1) + (0 - positions2D(2,t-1)) / (positions2D(2,t) - positions2D(2,t-1)) * dt;
                        positions2D(:,t) = positions2D(:,t-1) + (positions2D(:,t) - positions2D(:,t-1)) * (t_impact - time(t-1)) / dt;
                        positions2D(2,t) = 0;
                        hit_ground = true;
                    end
                else
                    positions2D(:,t) = positions2D(:,t-1);
                end
                
                % Sin air drag
                if ~hit_ground_no_drag
                    v_no_drag = [vx; vz - g * time(t)];
                    positions2D_no_drag(:,t) = positions2D_no_drag(:,t-1) + v_no_drag * dt;
                    
                    if positions2D_no_drag(2,t) <= 0
                        % Sacar el tiempo exacto de impacto
                        t_impact = time(t-1) + (0 - positions2D_no_drag(2,t-1)) / (positions2D_no_drag(2,t) - positions2D_no_drag(2,t-1)) * dt;
                        positions2D_no_drag(:,t) = positions2D_no_drag(:,t-1) + (positions2D_no_drag(:,t) - positions2D_no_drag(:,t-1)) * (t_impact - time(t-1)) / dt;
                        positions2D_no_drag(2,t) = 0;
                        final_x_no_drag = positions2D_no_drag(1,t);
                        hit_ground_no_drag = true;
                    end
                else
                    positions2D_no_drag(:,t) = [final_x_no_drag; 0];
                end
            end
            
            % Cortar arrays para tener solo datos válidos
            valid_points = max(find(positions2D(2,:) >= 0, 1, 'last'), find(positions2D_no_drag(2,:) >= 0, 1, 'last'));
            positions2D = positions2D(:, 1:valid_points);
            positions2D_no_drag = positions2D_no_drag(:, 1:valid_points);
            time_values = time(1:valid_points);
        end
        
        % Grafica los resultados de la simulación 2D
        % Muestra las trayectorias con y sin resistencia del aire
        function plot2DResults(app, positions2D, positions2D_no_drag, time_values)
            % Limpiar gráfico anterior
            cla(app.UIAxes2D);
        
            % Graficar trayectorias
            plot(app.UIAxes2D, positions2D(1,:), positions2D(2,:), 'LineWidth', 2, 'Color', 'r');
            hold(app.UIAxes2D, 'on');
            plot(app.UIAxes2D, positions2D_no_drag(1,:), positions2D_no_drag(2,:), 'LineWidth', 2, 'Color', 'b');
        
            % Etiquetas y título
            xlabel(app.UIAxes2D, 'Distance (m)');
            ylabel(app.UIAxes2D, 'Height (m)');
            title(app.UIAxes2D, 'Particle simulation');
            legend(app.UIAxes2D, 'with air drag', 'no air drag', 'Location', 'northeast');
        
            % Límites en ejes
            max_x = max(max(positions2D(1,:)), max(positions2D_no_drag(1,:)));
            max_y = max(max(positions2D(2,:)), max(positions2D_no_drag(2,:)));
            xlim(app.UIAxes2D, [0, max(100, max_x * 1.1)]);
            ylim(app.UIAxes2D, [0, max(50, max_y * 1.1)]);
        
            grid(app.UIAxes2D, 'on');
        
            % Datos de tabla
            table_data = [time_values', positions2D', positions2D_no_drag'];
            
            % Max. 50 elementos en tabla
            if size(table_data, 1) > 50
                indices = round(linspace(1, size(table_data, 1), 50));
                table_data = table_data(indices, :);
            end
            
            table_header = {'Time (s)', 'X (m)', 'Y (m)', 'X no drag (m)', 'Y no drag (m)'};
            uitable(app.UIFigure, 'Data', table_data, 'ColumnName', table_header, 'Position', [10 10 300 200], 'RowName', []);
            
            hold(app.UIAxes2D, 'off');
        end
        
        % Simulación 3D de la trayectoria de la partícula sobre el modelo del volcán
        % Calcula la trayectoria considerando la topografía del volcán
        function [positions3D, x, y, z] = run3DSimulation(app, dt, p, n, c, diameter, mass, shot_angle, initial_velocity)
            g = 9.81;
            total_time = 20; % Tiempo total de simulación (20s)
        
            % Crear modelo del volcán con funciones gaussianas
            [x, y] = meshgrid(-60000:1000:60000, -60000:1000:60000);
            z1 = 1117 * exp(-0.001*((x/1000+30).^2 + (y/1000+20).^2)); % Kitadake (Norte)
            z2 = 1060 * exp(-0.001*((x/1000-30).^2 + (y/1000-20).^2)); % Minamidake (Sur)
            z = z1 + z2;
        
            % Punto de erupción (Minamidake)
            eruption_x = 30000;
            eruption_y = 20000;
            eruption_z = interp2(x, y, z, eruption_x, eruption_y);
        
            % Convertir ángulo de tiro
            shot_angle_rad = deg2rad(shot_angle);
        
            % Componentes iniciales de velocidad (in m/s)
            vx = initial_velocity * cos(shot_angle_rad);
            vy = 0; % No hay velocidad en y
            vz = initial_velocity * sin(shot_angle_rad);
        
            % Calcular D/m (drag factor)
            A = pi * (diameter/2)^2;
            k = 0.5 * p * c * A / mass;
        
            % Arrays de simulación
            time = 0:dt:total_time;
            num_steps = length(time);
            positions3D = zeros(3, num_steps);
            velocities3D = zeros(3, num_steps);
            
            positions3D(:,1) = [eruption_x; eruption_y; eruption_z];
            velocities3D(:,1) = [vx; vy; vz];

            % Lógica de simulación igual a la 2D, con Método de Euler
            for t = 2:num_steps
                v = velocities3D(:,t-1);
                p = positions3D(:,t-1);
                speed = norm(v);
                
                % Calcular aceleración
                a_drag = -k * speed^(n-1) * v;
                a = a_drag + [0; 0; -g];
                
                % Actualizar velocidades y posiciones usando Método de Euler
                velocities3D(:,t) = v + a * dt;
                positions3D(:,t) = p + velocities3D(:,t) * dt;
            end
            
            % Cortar arrays a solo datos útiles
            positions3D = positions3D(:, 1:t);
        end
        
        % Grafica los resultados de la simulación 3D
        % Muestra la trayectoria sobre el modelo 3D del volcán
        function plot3DResults(app, positions, x, y, z)
            % Limpiar gráfico anterior
            cla(app.UIAxes3D);
            
            % Graficar superficie del volcán
            surf(app.UIAxes3D, x/1000, y/1000, z, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
            hold(app.UIAxes3D, 'on');
            colormap(app.UIAxes3D, 'hot');
            
            % Factor de escala para la trayectoria (solo componentes x e y)
            scale_factor = 100;
            
            % Calcular trayectoria escalada
            scaled_positions = positions;
            scaled_positions(1:2,:) = scaled_positions(1:2,:) - positions(1:2,1);  % Center the trajectory
            scaled_positions(1:2,:) = scaled_positions(1:2,:) * scale_factor;  % Scale only x and y
            scaled_positions(1:2,:) = scaled_positions(1:2,:) + positions(1:2,1);  % Move back to original position
            
            % Graficar trayectoria 3D escalada
            plot3(app.UIAxes3D, scaled_positions(1,:)/1000, scaled_positions(2,:)/1000, scaled_positions(3,:), 'LineWidth', 2, 'Color', 'r');
            
            % Etiquetas y títulos
            xlabel(app.UIAxes3D, 'X (km)');
            ylabel(app.UIAxes3D, 'Y (km)');
            zlabel(app.UIAxes3D, 'Z (m)');
            title(app.UIAxes3D, 'Particle simulation - 3D Sakurajima (xy-distance 100x)');
            
            % Aspect ratio y vista de gráfica 3D
            pbaspect(app.UIAxes3D, [1 1 0.5]);
            view(app.UIAxes3D, 45, 30);
            
            % Calcular límites
            max_height = max(max(positions(3,:)), max(z(:)));
            max_horizontal = max(max(abs(scaled_positions(1,:)-positions(1,1)), abs(scaled_positions(2,:)-positions(2,1)))) / 1000;
        
            % Límites de ejes para el gráfico 3D
            axis(app.UIAxes3D, 'tight');
            xlim(app.UIAxes3D, [min(x(:))/1000, max(x(:))/1000] + [-max_horizontal, max_horizontal]);
            ylim(app.UIAxes3D, [min(y(:))/1000, max(y(:))/1000] + [-max_horizontal, max_horizontal]);
            zlim(app.UIAxes3D, [0, max_height * 1.1]);
            
            % Poner cuadrícula/grid y hacer hold
            grid(app.UIAxes3D, 'on');
            hold(app.UIAxes3D, 'off');
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 1000 800];
            app.UIFigure.Name = 'Sakurajima particle simulation';

            % Create UIAxes3D
            app.UIAxes3D = uiaxes(app.UIFigure);
            app.UIAxes3D.Position = [320 60 660 520];

            % Create UIAxes2D
            app.UIAxes2D = uiaxes(app.UIFigure);
            app.UIAxes2D.Position = [320 600 660 180];

            % Create SimulateButton
            app.SimulateButton = uibutton(app.UIFigure, 'push');
            app.SimulateButton.ButtonPushedFcn = createCallbackFcn(app, @SimulateButtonPushed, true);
            app.SimulateButton.Position = [50 400 150 22];
            app.SimulateButton.Text = 'Simulate 2D & 3D';

            % Create TimestepEditFieldLabel
            app.TimestepEditFieldLabel = uilabel(app.UIFigure);
            app.TimestepEditFieldLabel.HorizontalAlignment = 'right';
            app.TimestepEditFieldLabel.Position = [0 741 120 22];
            app.TimestepEditFieldLabel.Text = 'Time step';

            % Create TimestepEditField
            app.TimestepEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimestepEditField.Position = [130 741 50 22];
            app.TimestepEditField.Value = 0.1;

            % Create pEditFieldLabel
            app.pEditFieldLabel = uilabel(app.UIFigure);
            app.pEditFieldLabel.HorizontalAlignment = 'right';
            app.pEditFieldLabel.Position = [71 700 25 22];
            app.pEditFieldLabel.Text = 'p';

            % Create pEditField
            app.pEditField = uieditfield(app.UIFigure, 'numeric');
            app.pEditField.Position = [111 700 50 22];
            app.pEditField.Value = 1.225;

            % Create nEditFieldLabel
            app.nEditFieldLabel = uilabel(app.UIFigure);
            app.nEditFieldLabel.HorizontalAlignment = 'right';
            app.nEditFieldLabel.Position = [71 659 25 22];
            app.nEditFieldLabel.Text = 'n';

            % Create nEditField
            app.nEditField = uieditfield(app.UIFigure, 'numeric');
            app.nEditField.Position = [111 659 50 22];
            app.nEditField.Value = 1;

            % Create cEditFieldLabel
            app.cEditFieldLabel = uilabel(app.UIFigure);
            app.cEditFieldLabel.HorizontalAlignment = 'right';
            app.cEditFieldLabel.Position = [71 618 25 22];
            app.cEditFieldLabel.Text = 'c';

            % Create cEditField
            app.cEditField = uieditfield(app.UIFigure, 'numeric');
            app.cEditField.Position = [111 618 50 22];
            app.cEditField.Value = 0.47;

            % Create DiameterofprojectileEditFieldLabel
            app.DiameterofprojectileEditFieldLabel = uilabel(app.UIFigure);
            app.DiameterofprojectileEditFieldLabel.HorizontalAlignment = 'right';
            app.DiameterofprojectileEditFieldLabel.Position = [0 577 120 22];
            app.DiameterofprojectileEditFieldLabel.Text = 'Diameter (m) ';

            % Create DiameterofprojectileEditField
            app.DiameterofprojectileEditField = uieditfield(app.UIFigure, 'numeric');
            app.DiameterofprojectileEditField.Position = [130 577 50 22];
            app.DiameterofprojectileEditField.Value = 0.6;

            % Create MassofprojectileEditFieldLabel
            app.MassofprojectileEditFieldLabel = uilabel(app.UIFigure);
            app.MassofprojectileEditFieldLabel.HorizontalAlignment = 'right';
            app.MassofprojectileEditFieldLabel.Position = [0 536 120 22];
            app.MassofprojectileEditFieldLabel.Text = 'Mass of proj. (kg)';

            % Create MassofprojectileEditField
            app.MassofprojectileEditField = uieditfield(app.UIFigure, 'numeric');
            app.MassofprojectileEditField.Position = [130 536 50 22];
            app.MassofprojectileEditField.Value = 1;

            % Create ShotangleEditFieldLabel
            app.ShotangleEditFieldLabel = uilabel(app.UIFigure);
            app.ShotangleEditFieldLabel.HorizontalAlignment = 'right';
            app.ShotangleEditFieldLabel.Position = [0 495 120 22];
            app.ShotangleEditFieldLabel.Text = 'Shot angle (°)';

            % Create ShotangleEditField
            app.ShotangleEditField = uieditfield(app.UIFigure, 'numeric');
            app.ShotangleEditField.Position = [130 495 50 22];
            app.ShotangleEditField.Value = 40;

            % Create InitialvelocityEditFieldLabel
            app.InitialvelocityEditFieldLabel = uilabel(app.UIFigure);
            app.InitialvelocityEditFieldLabel.HorizontalAlignment = 'right';
            app.InitialvelocityEditFieldLabel.Position = [0 454 120 22];
            app.InitialvelocityEditFieldLabel.Text = 'Initial velocity (m/s)';

            % Create InitialvelocityEditField
            app.InitialvelocityEditField = uieditfield(app.UIFigure, 'numeric');
            app.InitialvelocityEditField.Position = [130 454 50 22];
            app.InitialvelocityEditField.Value = 50;
        end
    end

    methods (Access = public)

        % Construct app
        function app = SakurajimaApp
            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end