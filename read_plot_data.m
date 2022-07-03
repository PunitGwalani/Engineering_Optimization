format long;

% delta_v=dlmread("Results101/Final_DeltaV_og.dat");
tof=dlmread("Results101/All_files/AGAIN_TimeOfFlight.dat");
epoch=dlmread("Results101/All_files/AGAIN_DepartureEpochs.dat");
delta_v2=dlmread("Results101/All_files/AGAIN_Venus_Capture_DeltaV.dat");


%Gravitational acceleration on Earth (m/s2)
g_E=9.81;

%AU in m
au = 1.495978 * 10^(11);
max_limit = 1.25; %AU
min_limit = 0.56; %AU

max_pos = (dlmread("Results101/All_files/AGAIN_Max_positions.dat"))/au;
min_pos = (dlmread("Results101/All_files/AGAIN_Min_positions.dat"))/au;

%Initial mass of propellant + dry mass (kg)
M_init=1245; 
M_fuel = 570;

%Specific Impulse of rocket engine (s) [Currently based on Venus Express]
I_sp=317;

%km/s to m/s
% delta_v=delta_v'*1000;
% delta_v1=delta_v1*1000;
% delta_v2=delta_v2*1000;

constraint_delta_v = (I_sp*g_E*log(M_init/M_fuel))/1000

activate_mass_penalties =1;
activate_max_pos_penalties = 1;
activate_min_pos_penalties = 1;

x = [];
y = [];

for i=1:300
    for j=1:150
        
        if activate_mass_penalties == 1
            if delta_v2(i,j) > constraint_delta_v
                x = [x; epoch(i,j)];
                y = [y; tof(i,j)];
%                 delta_v2(i,j) = delta_v2(i,j) + 1e30;
            end
        end
        
        if activate_max_pos_penalties == 1
            if max_pos(i,j) > max_limit
                x = [x; epoch(i,j)];
                y = [y; tof(i,j)];
%                 delta_v2(i,j) = delta_v2(i,j) + 10000*abs(max_pos(i,j) - max_limit);
            end
        end
        
        if activate_min_pos_penalties == 1
            if min_pos(i,j) < min_limit
                x = [x; epoch(i,j)];
                y = [y; tof(i,j)];
%                 delta_v2(i,j) = delta_v2(i,j) + 10000*abs(min_pos(i,j) - min_limit);
            end
        end
                   
    end
end

epoch_julian = linspace(min(epoch(1,:)), max(epoch(1,:)), 7);
epoch_datetime = string(datetime(epoch_julian + 2451545.0,'convertfrom','juliandate','Format','yyy-MMM-dd'));

[C,h] = contourf(epoch(1,:),tof(:,1),delta_v2,...
    [0.03, 0.05,0.1 ,0.5, 1, 2, 4, 10, 30, 50, 60, 500],...
'ShowText','on')
% colormap('turbo')
clabel(C,h,'FontSize',9,'Color','white');

xticks(epoch_julian)
xticklabels(epoch_datetime)
zlabel('\Delta V [km/s]')
a = colorbar;
a.Label.String = '\Delta V [km/s]';
xlabel('Launch Dates [2024]','FontWeight','bold')
ylabel('Time of Flight [Days]','FontWeight','bold')

hold on
sz = 5;
scatter(x,y,sz,'MarkerFaceColor','r', 'MarkerEdgeColor', '#D95319');

% Contour_edits(epoch(1,:),tof(:,1),delta_v2)
% surfc(epoch(1,:),tof(:,1),delta_v2)

% 
% hYLabel = get(gca,'YLabel');
% set(hYLabel,'rotation',-30,'VerticalAlignment','middle')
% 
% hXLabel = get(gca,'XLabel');
% set(hXLabel,'rotation',20,'HorizontalAlignment','center')


% contour(epoch(1,:),tof(:,1),delta_v2/1000,20,'ShowText','on')


