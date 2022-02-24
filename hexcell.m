function [D] = hexcell(SP)

Nu = SP.Nu; %Number of users
Nc = SP.Nc; 
dmax = SP.dmax; %Radius of Hexagon
dmin = SP.dmin;
d0 = SP.d0; % meter
rho = SP.rho; % pathloss exponent
s = SP.s; % (dB) shadow fading standard deviation
LD = SP.LD;
NF = SP.NF;
BW = SP.BW;


CX = sqrt(3)*dmax * [0, cos((0:Nc-2)*pi/3+pi/6)]; % x coordinates of centers
CY = sqrt(3)*dmax * [0, sin((0:Nc-2)*pi/3+pi/6)]; % y coordinates of centers



points = zeros(Nu,Nc,2);


    for c = 1:Nc
        
        v_x_max = dmax * cos((0:6)*pi/3) + CX(c); % vertax of hexagon centered at center
        v_y_max = dmax * sin((0:6)*pi/3) + CY(c);

        v_x_min = dmin * cos((0:6)*pi/3) + CX(c);
        v_y_min = dmin * sin((0:6)*pi/3) + CY(c);
        %% Place Nu samples in the hexagon
        
        [p_x, p_y] = makepoints(SP);
        
        x = CX(c) + p_x;
        y = CY(c) + p_y;
        
        %% Coordinates of Nu points in a cell
        
        points(:,c,1) = x;
        points(:,c,2) = y;
        
%         hold on
%         plot(CX(c), CY(c), 'o');
%         plot(v_x_max, v_y_max);
%         plot(v_x_min, v_y_min);
%         plot(x, y, '.');
        
    end
    
    D = zeros(Nu,Nc,Nc);
    
    for c = 1:Nc
       diff_x = points(:,:,1) - [CX(c)]; % distance is measured between corresponding cell center
       diff_y = points(:,:,2) - [CY(c)];

       D(:,:,c) = sqrt(diff_x.^2 + diff_y.^2); % D(:,:,c)
    end
    
end

 




