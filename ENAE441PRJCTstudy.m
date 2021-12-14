% Jacob Dickinson-Sabonis
% Implenting code written by Alexander 
clear
lat = -30.142803; %deg
long = -70.694528; %deg
alt = 1500 ; %m
mu = 3.986e5;
re = 6.378e3;
load('opt2satCset3');
load('opt3satCset3');
set1 = opt2satCset3([50 70 90],:); %bestest
set2 = opt2satCset3([120 234 567],:);
set3 = opt2satCset3([20 80 150], :);
set4 = opt2satCset3([450 480 520],:);
set5 = opt2satCset3([800 840 880],:); %best
set6 = opt2satCset3([760 830 900],:);
set7 = opt2satCset3([600 630 660],:);

obs1 = [1:5:100]';
obs2 = [101:25:600]';
obs3 = [1:10:200]';
obs4 = [401:10:600]';
obs5 = [801:5:900]';
obs6 = [701:10:900]';
obs7 = [601:5:700]';
observations = [obs1 obs2 obs3 obs4 obs5 obs6 obs7];

%from Gauss --> Gibbs
posvel1 = [39872038.9674329;11590157.3549449;-125207.951618490;
 -869.578048678686;2995.80289919241;0.999402099155789];
posvel2 = [33956576.2110714;23126642.9144956;-155635.822122117;
    -1782.11675578726;2604.37443529279;1.87224698835830];
posvel3 = [39212399.8236516;13012258.2446220;-142537.808515821;
    -984.225578887384;2974.71608008337;1.34648278919514];
posvel4 = [21485268.0398972;35715876.3749781;-85499.1248084330; 
    -2668.57483724116;1602.01635853691;3.84337377253836];
posvel5 = [-9028546.27809327;40702422.0492338;-34590.4258942042; 
    -3032.13491491756;-672.880809412959;5.15846241104134];
posvel6 = [-8019568.92059573;40860283.9240002;-40725.8545175755; 
-3059.37502627744;-601.433284881025;5.18030678831449];
posvel7 = [9722717.57013852;40531203.9550206;-66772.3988078783;
    -3020.45624181846;725.192917058197;4.90145821474636];
posvels = [posvel1 posvel2 posvel3 posvel4 posvel5 posvel6 posvel7];


for z=1:length(posvels)
%% Varying Force Models
%posvel2 is the 2nd data set, now we are examining, input is 1,2,3
%corresponding to 2x0, 2x2, and 20x20 respectively 
%set the set 2 data, 
chile.lla = latlonalt_deg(lat,long,alt);
posvel = posvels(:,z);
observe = observations(:,z);
mu = 3.986e5;
re = 6.378e3;
load('opt2satCset3');
load('opt3satCset3');

% initial orbit determination (2 and 6)
set = opt2satCset3([50 70 90],:);
epoch = set.datetime(2);
state_init = pvt(epoch,posvel(1:3), posvel(4:6));
%changing solar radiation pressure 

force_model_2x0_1m2_cr1 = force_model(2,0,0,0,40,1.2,4000); % varied force model, little noticeable RMS change, 2x0 / 2x2 were lowest
oapchile = make_station("OAP-Chile", lat, long, alt);
fourcols = ["observation_number" "datetime"...
"azimuth_deg" "elevation_deg"];
night1 = opt2satCset3(observe,fourcols); % varied the observations (early/middle/late and number of points) 
od_night1 = determine_orbit(state_init, oapchile,...
night1, force_model_2x0_1m2_cr1);
sat_force.estimated = od_night1.estimated;
RMS_force = od_night1.details;
%input x = 2 for 2x0, no drag solar pressure model --> res = .3546
%% Orbit Validation of Estimation
clear resel resaz res
% validation using data from the 2nd night
smallset = opt3satCset3(10:40:1000,:);
force_model_2x0_1m2_val = force_model_third_body(20,20,0,0,4,1.2,4000);
sat_force.ephemeris = propagate_to_times(sat_force.estimated, smallset.datetime, force_model_2x0_1m2_val);
N = height(smallset);
sum = 0;
for i = 1:N
sat_force.aer(i) = eci_to_azelrn(sat_force.ephemeris.epoch(i),sat_force.ephemeris.position_m(i,:),chile.lla);
resaz(i) = smallset.azimuth_deg(i) - sat_force.aer(i).azimuth_deg;
resel(i) = smallset.elevation_deg(i) - sat_force.aer(i).elevation_deg;
sum = sum + resaz(i)^2+resel(i)^2;
end
res_force(z) = sqrt(1/N*sum)
dtepoch = seconds(smallset.datetime-smallset.datetime(1));
%{
figure(1)
plot(dtepoch,resaz)
hold on
title('az residual (varying force)')
figure(2)
plot(dtepoch, resel)
hold on 
title('el residual (varying force)')
%}

%% Varying Data while keeping force model constant 
% checks RMS while keeping the force model constant at 2x0

set = set2; %varied the input to determine optimal data set
epoch = set.datetime(2);
state_init = pvt(epoch,posvel(1:3), posvel(4:6));
% Orbit Estimation
force_model_2x0 = force_model_third_body(2,0,0,0,40,1.2,4000);
oapchile = make_station("OAP-Chile", lat, long, alt);
fourcols = ["observation_number" "datetime"...
"azimuth_deg" "elevation_deg"]; 
night1_early_25pts = opt2satCset3(observe,fourcols);
od_night1_early = determine_orbit(state_init, oapchile,...
night1_early_25pts, force_model_2x0);
sat_data.estimated = od_night1_early.estimated;
RMS = od_night1_early.details;

%% Orbit Validation of Estimation while varying data
clear resel resaz res
% validation using data from the 2nd night
smallset = opt3satCset3(10:40:1000,:);
force_model_2x0_1m2_val = force_model_third_body(20,20,0,0,4,1.2,4000);
sat_data.ephemeris = propagate_to_times(sat_data.estimated, smallset.datetime, force_model_2x0_1m2_val);
N = height(smallset);
sum = 0;
sum = 0;
for i = 1:N
sat_data.aer(i) = eci_to_azelrn(sat_data.ephemeris.epoch(i),sat_data.ephemeris.position_m(i,:),chile.lla);
resaz(i) = smallset.azimuth_deg(i) - sat_data.aer(i).azimuth_deg;
resel(i) = smallset.elevation_deg(i) - sat_data.aer(i).elevation_deg;
sum = sum + resaz(i)^2+resel(i)^2;
end
res_data(z) = sqrt(1/N*sum)

%{
dtepoch = seconds(smallset.datetime-smallset.datetime(1));
figure(1)
plot(dtepoch,
az)
hold on
title('az residual (varying data)')
figure(2)
plot(dtepoch, resel)
hold on 
title('el residual (varying data)')
%}
end
