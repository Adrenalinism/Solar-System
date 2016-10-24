clear;
clc;
clf;
diary on
hold on

earth_moon_radius_au = 1.16138017*10^(-5);
tilt = 0.4;
largest_orbit = 39.5;

% System Type
if randi([0 1],1,1)  
    system_type = 'Binary';
else
    system_type = 'Singular';
end

% Star(s) Mass (relative to Sun)
star_mass = gevrnd(0.5,0.6,1.1);
if star_mass < 0.55
    star_mass = 0.55;
end
if strcmp(system_type,'Binary')
    star_2_mass = gevrnd(0.5,0.6,1.1);
    if star_2_mass < 0.55
        star_2_mass = 0.55;
    end
end
if strcmp(system_type,'Binary')
    correlated_star_mass = 1.5*((star_mass+star_2_mass)/2);
else
    correlated_star_mass = star_mass;
end

% Habitable Zone
habitable_radius_inner = sqrt(correlated_star_mass^3.5/1.1);
habitable_radius_outer = sqrt(correlated_star_mass^3.5/0.53);
if habitable_radius_outer > largest_orbit
    largest_orbit = habitable_radius_outer;
end
plot(habitable_radius_inner*cos(0:0.01:2*pi),habitable_radius_inner*tilt*sin(0:0.01:2*pi),'--w');
plot(habitable_radius_outer*cos(0:0.01:2*pi),habitable_radius_outer*tilt*sin(0:0.01:2*pi),'--w');
plot(39.5*cos(0:0.01:2*pi),39.5*tilt*sin(0:0.01:2*pi),':w');

% Number of Planets in System
num_planets = round((0.3*correlated_star_mass+0.7)*normrnd(6,2.5));
if num_planets <= 0
    num_planets = 1;
end

% CMD Displays
disp('GENERATING SOLAR SYSTEM.....')
disp(' ');
disp(['System_Type: ',system_type]);
if strcmp(system_type,'Binary')
    disp(['Star Masses: ' num2str(star_mass), ' and ', num2str(star_2_mass) ' relative to Sun ']);
else
    disp(['Star Mass: ', num2str(star_mass),' relative to Sun ']);
end
disp(['Number of planets: ', num2str(num_planets)]);
disp(' ');
disp(' ');

%-----------------------------PLANETS--------------------------------------
planet_counter = 0;
while (planet_counter < num_planets)    

    % Planet Type
    if randi([0 9],1,1) <= 5    
        planet_type = 'Rocky';
    else
        planet_type = 'Gassy';
    end
    disp(['Planet ',num2str(planet_counter+1)]);
    disp(['Planet Type: ',planet_type]);

    % Planet Mass (relative to Earth)
    if strcmp(planet_type,'Gassy')
        planet_mass = randi([8 450],1,1);
    else
        planet_mass = normrnd(1.2,0.3);
        if planet_mass <= 0.03
            planet_mass = 0.04675;
        end
    end
    disp(['Planet Size: x', num2str(planet_mass), ' of earth mass']);
    

    % Planet Distance (in AU)
    planet_distance_trizone = randi([0 2],1,1);
    if planet_distance_trizone == 0
        planet_distance_au = (rand + 0.5)*(0.3*correlated_star_mass+0.7);
    elseif planet_distance_trizone == 1
        planet_distance_au = ((randi([15000 100000],1,1))/10000)*(0.3*correlated_star_mass+0.7);
    else
        planet_distance_au = ((randi([100000 450000],1,1))/10000)*(0.3*correlated_star_mass+0.7);    
    end
    if planet_distance_au > largest_orbit
        largest_orbit = planet_distance_au;
    end
    disp(['Planet Distance: ', num2str(planet_distance_au), ' AU']);
    
    % Planet Period
    planet_distance_meters = 1.496e+11*planet_distance_au;
    correlated_star_mass_kg = 1.989*10^30*correlated_star_mass;
    planet_period = (sqrt((planet_distance_meters^3*4*pi*pi)/(6.67*10^(-11)*correlated_star_mass_kg))/(60*60*8760));
    disp(['Planetary Period: ', num2str(planet_period), ' earth years']);

    % Planet Number of Moons
    if planet_mass <=  4    
       moon_quadzone = randi([1 1000],1,1);
       if moon_quadzone <= 500
           planet_moons = 0;
       elseif moon_quadzone > 500 && moon_quadzone <= 850
           planet_moons = 1;
       else
           planet_moons = 2;
       end
    elseif planet_mass > 4 && planet_mass <= 100
        planet_moons = round((round(0.15*planet_mass + 2.55) + randi([0 4],1,1) - 2)/8);
    elseif planet_mass > 100 && planet_mass <= 225
        planet_moons = round((round(0.15*planet_mass + 2.55) + randi([0 18],1,1) - 9)/8);
    else
        planet_moons = round((round(0.15*planet_mass + 2.55) + randi([0 40],1,1) - 20)/8);    
    end
    disp(['Number of Major Moons: ' num2str(planet_moons)]);

    % Surface Water
    water_source_chance = randi([0 1],1,1); 
    surface_water_gate = 0;
    too_hot = '';
    too_cold = '';
    no_solid_ground = '';
    no_water_source = '';
    
    if planet_distance_au <= habitable_radius_inner
        too_hot = 'too hot, ';
        surface_water_gate = 1;
    end
    if planet_distance_au >= habitable_radius_outer
        too_cold = 'too cold, ';
        surface_water_gate = 1;
    end
    if strcmp(planet_type,'Gassy')  
        no_solid_ground = 'no solid ground, ';
        surface_water_gate = 1;
    end
    if water_source_chance == 1
        no_water_source = 'no H2O source';
        surface_water_gate = 1;
    end
    if surface_water_gate == 1
        if water_source_chance == 0
            surface_water = (['Surface Water: False (',too_hot,too_cold,no_solid_ground,no_water_source, ')']);
            surface_water = surface_water(1:end-3);
            surface_water = strcat(surface_water,')');
            disp(surface_water);
        else
            disp(['Surface Water: False (',too_hot,too_cold,no_solid_ground,no_water_source, ')']);
        end
    end

    if planet_distance_au >= habitable_radius_inner && planet_distance_au <= habitable_radius_outer && strcmp(planet_type,'Rocky') && water_source_chance == 0
        disp('Surface Water: True (GOLDILOCKS!!)');
    end
    
    % Planet Radius + Plot Planet Body
    if strcmp(planet_type,'Gassy')
        planet_density = normrnd(1.3,0.15);
    else
        planet_density = normrnd(5,0.3);
    end
    planet_randomizer = 2*pi*rand(1,1);
    planet_radius_au = (((nthroot((3*planet_mass*5.972*10^24)/(planet_density*4*pi), 3)))/1000)*6.68459e-9;
    rectangle('Position',[(planet_distance_au*cos(planet_randomizer)-planet_radius_au) (planet_distance_au*tilt*sin(planet_randomizer)-planet_radius_au) (2*planet_radius_au) (2*planet_radius_au)],'Curvature',[1 1],'EdgeColor','w','FaceColor','w');
    
    % Plot Planet Orbit
    plot(planet_distance_au*cos(0:0.01:2*pi),planet_distance_au*tilt*sin(0:0.01:2*pi));
    if planet_distance_au >= habitable_radius_inner && planet_distance_au <= habitable_radius_outer && strcmp(planet_type,'Rocky') && water_source_chance == 0
        plot(planet_distance_au*cos(planet_randomizer),planet_distance_au*tilt*sin(planet_randomizer),'xc');
    else
        plot(planet_distance_au*cos(planet_randomizer),planet_distance_au*tilt*sin(planet_randomizer),'xw');
    end   
    disp(' ')
    planet_counter=planet_counter+1;
    
    % Plot Moons
    moon_counter = 0;
    while moon_counter < planet_moons       
        moon_radius_au = normrnd(0.8*earth_moon_radius_au,0.25*earth_moon_radius_au);
        if moon_radius_au < earth_moon_radius_au*0.1
            moon_radius_au = earth_moon_radius_au*0.1;
        end
        moon_distance_au = ((randi([150 3000],1,1))/100)*planet_radius_au;
        rectangle('Position',[(planet_distance_au*cos(planet_randomizer)-moon_distance_au-moon_radius_au) (planet_distance_au*tilt*sin(planet_randomizer)-moon_radius_au) (2*moon_radius_au) (2*moon_radius_au)],'Curvature',[1 1],'EdgeColor','w','FaceColor','w');
        
        plot(moon_distance_au*cos(0:0.01:2*pi)+(planet_distance_au*cos(planet_randomizer)),moon_distance_au*tilt*sin(0:0.01:2*pi)+(planet_distance_au*tilt*sin(planet_randomizer)),'Color',[0.16 0.16 0.16],'LineWidth',0.02);

        moon_counter = moon_counter + 1;
    end
end
%--------------------------------------------------------------------------

% Plot Star(s)
if strcmp(system_type,'Singular')
    star_radius = (-0.0123*correlated_star_mass*correlated_star_mass + 0.6319*correlated_star_mass + 0.4829)*0.00464913034;
    rectangle('Position',[-star_radius -star_radius 2*star_radius 2*star_radius],'Curvature',[1 1],'EdgeColor','w','FaceColor','w');
else
    star_radius = (-0.0123*star_mass*star_mass + 0.6319*star_mass + 0.4829)*0.00464913034;
    star_radius_2 = (-0.0123*star_2_mass*star_2_mass + 0.6319*star_2_mass + 0.4829)*0.00464913034;
    
    binary_star_locations = randi([0 2],1,1) - 1;
    while binary_star_locations == 0
        binary_star_locations = randi([0 2],1,1) - 1;
    end   
    star_offset = (randi([150 280],1,1)/10000);
    rectangle('Position',[(-star_radius-star_offset)*binary_star_locations -star_radius+star_offset 2*star_radius 2*star_radius],'Curvature',[1 1],'EdgeColor','w','FaceColor','w');
    rectangle('Position',[(-star_radius_2+star_offset)*binary_star_locations -star_radius_2-star_offset 2*star_radius_2 2*star_radius_2],'Curvature',[1 1],'EdgeColor','w','FaceColor','w');
end
plot(0.00464913*cos(0:0.01:2*pi),0.00464913*sin(0:0.01:2*pi),'m');

% Figure Colors/Labelling and CMD outputs
set(gca,'Color','k')
set(gca,'xcolor','k') 
set(gca,'ycolor','k') 
% system_name = strcat('Star System: ',strcat(system_type,strcat(' Sun Mass Star, ',strcat(strcat(num2str(planet_counter))))));
% title(system_name)
xlabel('Distance in AU (0,0 is centre of system)')
axis([-(2+largest_orbit) (2+largest_orbit) -(2+largest_orbit) (2+largest_orbit)])
disp(' ')
disp('=======================================================================================================')
disp(' ')
diary off
hold off


    




