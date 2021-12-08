clear
mu = 3.986e5;
re = 6.378e3;
%% Load Sat data
load('opt2satCset3');
load('opt3satCset3');

%% Initial Orbit Determination

set1 = opt2satCset3([5 15 30],:);
set2 = opt2satCset3([50 70 90],:);
set3 = opt2satCset3([120 140 160],:);
set4 = opt2satCset3([100 400 700], :);
set5 = opt2satCset3([200 525 850],:);
set6 = opt2satCset3([850 870 890],:); 
set7 = opt2satCset3([885 900 915],:);

sets = [set1; set2; set3; set4; set5; set6; set7];
numsets = size(sets,1);

lat = sets.site_latitude_deg(1);
long = sets.site_longitude_deg(1);
alt = set1.site_altitude_m(1);
chile.lla = latlonalt_deg(lat, long, alt);
for i = 0:numsets/3-1
    % Gibbs
    R = zeros(3,3);
    L = zeros(3,3);
    t2 = sets.datetime(3*i+2);
    for j = 1:3
        index = 3*i+j;
        chile.lla.epoch = sets.datetime(index);
        chile.eci = eci(chile.lla);
        R(:,j) = chile.eci.position_m'*1e-3;
        raant = sets.right_ascension_deg(index);
        dect = sets.declination_deg(index);
        L(:,j) = [cosd(dect)*cosd(raant); cosd(dect)*sind(raant); sind(dect)];
        tau(j) = seconds(chile.lla.epoch - t2);
    end
    M = L\R;
    a1 = tau(3)/(tau(3)-tau(1));
    a3 = -tau(1)/(tau(3)-tau(1));
    a1u = tau(3)*((tau(3)-tau(1))^2-tau(3)^2)/(6*(tau(3)-tau(1)));
    a3u = -tau(1)*((tau(3)-tau(1))^2-tau(1)^2)/(6*(tau(3)-tau(1)));
    A = M(2,1)*a1-M(2,2)+M(2,3)*a3;
    B = M(2,1)*a1u+M(2,3)*a3u;
    E = dot(L(:,2),R(:,2));
    R2_2 = norm(R(:,2))^2;
    poly = [1 0 -(A^2+2*A*E+R2_2) 0 0 2*mu*B*(A+E) 0 0 -mu^2*B^2];
    r2_1 = roots(poly);
    r2_1 = r2_1(imag(r2_1)==0);
    r2_1 = r2_1(r2_1>re);
    u = mu/r2_1^3;
    c1 = -(-a1-a1u*u);
    c2 = -1;
    c3 = -(-a3-a3u*u);
    C = [-c1 -c2 -c3]';
    Rho = (M*C)./(-C);
    r1 = R(:,1) + Rho(1)*L(:,1);
    r2 = R(:,2) + Rho(2)*L(:,2);
    r3 = R(:,3) + Rho(3)*L(:,3);
    % Gibbs
    D = cross(r2,r3)+cross(r3,r1)+cross(r1,r2);
    N = cross(norm(r1)*r2,r3)+cross(norm(r2)*r3,r1) + cross(norm(r3)*r1,r2);
    S = (norm(r2)-norm(r3))*r1 + (norm(r3) - norm(r1))*r2 + (norm(r1)-norm(r2))*r3;
    W = N/norm(N);
    Q = S/norm(S);
    P = cross(Q,W);
    e = norm(S)/norm(D)
    p = norm(N)/norm(D);
    a = p/(1-e^2)
    n = cross([0;0;1],W)/norm(cross([0;0;1],W));
    r2dot = sqrt(mu/(norm(N)*norm(D)))*(cross(D,r2/norm(r2))+S)
    oe(i+1,:) = rv2oe1(r2,r2dot,mu)
    posvel(i+1) = pvt(t2,r2*1e3,r2dot*1e3);
end

%% Validation of initial orbit determination
smallset = opt3satCset3(50:40:1000,:);
force_model_4x4_1m2_cr1 = force_model_third_body(4,4,0,0,1,1,1000);
for i = 1:numsets/3
    sat.ephemeris(i) = propagate_to_times(posvel(i), smallset.datetime, force_model_4x4_1m2_cr1);
end
%%
N = height(smallset);
for j = 1:numsets/3
    sum = 0;
    for i = 1:N
    sat.aer(i) = eci_to_azelrn(sat.ephemeris(j).epoch(i),sat.ephemeris(j).position_m(i,:),chile.lla);
    resaz(i,j) = smallset.azimuth_deg(i) - sat.aer(i).azimuth_deg;
    resel(i,j) = smallset.elevation_deg(i) - sat.aer(i).elevation_deg;
    sum = sum + resaz(i,j)^2+resel(i,j)^2;
    end
res(j) = sqrt(1/N*sum);
end
% Second Data Set Lowest RMS
%% Orbit Estimation
force_model_4x4_1m2_cr1 = force_model_third_body(4,4,0,0,1,1,1000);
oapchile = make_station("OAP-Chile", lat, long, alt);
fourcols = ["observation_number" "datetime"...
"azimuth_deg" "elevation_deg"];
night1_early_25pts = opt2satCset3([1:4:100],fourcols);
od_night1_early = determine_orbit(posvel(2), oapchile,...
night1_early_25pts, force_model_4x4_1m2_cr1)
sat2.estimated = od_night1_early.estimated;
%%
smallset = opt3satCset3(50:40:1000,:);
force_model_4x4_1m2_cr1 = force_model_third_body(4,4,0,0,1,1,1000);
sat2.ephemeris(i) = propagate_to_times(sat2.estimated, smallset.datetime, force_model_4x4_1m2_cr1);
%%
N = height(smallset);
for j = 1
    sum = 0;
    for i = 1:N
    sat2.aer(i) = eci_to_azelrn(sat2.ephemeris.epoch(i),sat2.ephemeris(j).position_m(i,:),chile.lla);
    resaz(i,j) = smallset.azimuth_deg(i) - sat2.aer(i).azimuth_deg;
    resel(i,j) = smallset.elevation_deg(i) - sat2.aer(i).elevation_deg;
    sum = sum + resaz(i,j)^2+resel(i,j)^2;
    end
res(j) = sqrt(1/N*sum);
end