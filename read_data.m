format long;

delta_v=dlmread("Final_DeltaV_og.dat");
tof=dlmread("Final_TimeOfFlight_og.dat");
epoch=dlmread("Final_DepartureEpochs_og.dat");
delta_v1=dlmread("Earth_Escape_DeltaV_final.dat");
delta_v2=dlmread("Venus_Capture_DeltaV_final.dat");

%Gravitational acceleration on Earth (m/s2)
g_E=9.81;

%Initial mass of propellant + dry mass (kg)
M_init=1245; 

%Specific Impulse of rocket engine (s) [Currently based on Venus Express]
I_sp=317;

%km/s to m/s
delta_v=delta_v'*1000;
delta_v1=delta_v1*1000;
delta_v2=delta_v2*1000;


for i=1:300
    for j=1:150

        %Mass ratio: Mass used during Venus Capture/Initial mass
        Mass_ratio(i,j)=exp(-(delta_v2(i,j)/(I_sp*g_E)))*(exp((delta_v2(i,j)/(I_sp*g_E)))-1);

        %Mass of propellant used during Earth escape
        %M_E(i,j)=exp(-(delta_v1(i,j)/(I_sp*g_E)))*M_init*(exp((delta_v1(i,j)/(I_sp*g_E)))-1);

        %Mass of propellant used during Venus capture
        M_V(i,j)=M_init*Mass_ratio(i,j);
        
    end
end

% surf(epoch(1,:),tof(:,1),delta_v)
% xlabel('Epoch')
% ylabel('Time of Flight')
contour(epoch(1,:),tof(:,1),delta_v/1000,[4 5 6 7 8 9 10 30 50 60 70 80 90 95],'ShowText','on')
%contour(epoch(1,:),tof(:,1),M_V,[0 800 1000 1200],'ShowText','on')


