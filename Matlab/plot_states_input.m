function [y] = plot_states_input(x_evol, tx, u_opt, tu)
    n = size(x_evol, 2);
    m = size(u_opt, 1);

    if n == 2
        if m == 2
            figure
            subplot(2, 2, 1)
            plot(tx, x_evol(:, 1))
            grid
            title('State 1')

            subplot(2, 2, 2)
            plot(tx, x_evol(:, 2))
            grid
            title('State 2')

            subplot(2, 2, 3)
            stem(tu, u_opt(1,:))
            ylim([-2, 2])
            grid
            title('Input 1')

            subplot(2, 2, 4)
            stem(tu, u_opt(2,:))
            ylim([-2, 2])
            grid
            title('Input 2')

            return
        end
    end

    if n == 4
        if m == 4
            figure
            subplot(2, 4, 1)
            plot(tx, x_evol(:, 1))
            grid
            title('State 1')
            
            subplot(2, 4, 2)
            plot(tx, x_evol(:, 2))
            grid
            title('State 2')
            
            subplot(2, 4, 3)
            plot(tx, x_evol(:, 3))
            grid
            title('State 3')
            
            subplot(2, 4, 4)
            plot(tx, x_evol(:, 4))
            grid
            title('State 4')
            
            subplot(2, 4, 5)
            stem(tu, u_opt(1,:))
            ylim([-2, 2])
            grid
            title('Input 1')
            
            subplot(2, 4, 6)
            stem(tu, u_opt(2,:))
            ylim([-2, 2])
            grid
            title('Input 2')
            
            subplot(2, 4, 7)
            stem(tu, u_opt(3,:))
            ylim([-2, 2])
            grid
            title('Input 1')
            
            subplot(2, 4, 8)
            stem(tu, u_opt(4,:))
            ylim([-2, 2])
            grid
            title('Input 2')
        end
    end

    if n == 5
        if m == 2
            figure
            subplot(2, 1, 1)
            plot(tx, x_evol(:, 1))
            grid
            title('Roll rate')

            subplot(2, 1, 2)
            plot(tx, x_evol(:, 2))
            grid
            title('Yaw rate')

            figure
            subplot(2, 1, 1)
            plot(tx, x_evol(:, 3))
            grid
            title('Roll angle')

            subplot(2, 1, 2)
            plot(tx, x_evol(:, 4))
            grid
            title('Yaw angle')

            figure
            plot(tx, x_evol(:, 5))
            grid
            title('Y Axis speed (body RS)')
        
            figure
            subplot(2, 1, 1)
            stem(tu, u_opt(1, :))
            grid
            title('Aileron angle (roll)')

            subplot(2, 1, 2)
            stem(tx, u_opt(2, :))
            grid
            title('Rudder angle (yaw)')
        end
    end

    if n == 2
        if m == 1
            figure
            subplot(2, 1, 1)
            plot(tx, x_evol(:, 1))
            grid
            title('Speed')
            xlabel('Seconds') 
            ylabel('Meters / second') 

            subplot(2, 1, 2)
            plot(tx, x_evol(:, 2))
            grid
            title('Position')
            xlabel('Seconds') 
            ylabel('Meters') 
            
            fontsize(16,"points")

            figure
            stem(tu, u_opt)
            grid
            title('Thrust')
            xlabel('Seconds') 
            ylabel('Thruster capacity (0% - 100%)') 

            fontsize(16,"points")
        end
    end
end


