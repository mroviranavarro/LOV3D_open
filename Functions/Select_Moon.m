function[Moon_cst] = Select_Moon(moon_ID)

%               ----------------------------
%                   FUNCTION Select_Moon
%               ----------------------------
% 
% INPUT:    - moon_ID:  number in which the first digit refers to the
%                       planet and the second to the satellite
%
% OUTPUT:   - Moon_cst: Structure containing information about the
%                       satellite's orbital parameters, period, radius and
%                       mass of the parent planet.
%
% Important: Function is in development and still needs data for several
% moons.
% Moon data can also be subject to updates according to more recent
% observations (Parameters listed below have been obtained from
% NASA/JPL's Solar System Dynamics database and [Seidelman et al., 2007]).

Moon_cst.G      = 6.67428*10^(-11);
Moon_cst.kyear  = 3.156*10^10;
Moon_cst.mum=170e9/(2*(1+0.33));
switch moon_ID
    % Earth's moon
    case 31
        % Moon    
        Moon_cst.Mp     = 59.736E23; 
        Moon_cst.r      = 1737.4E3; 
        
        Moon_cst.Mass=  7.3458e+22;
        Moon_cst.a      = 384400000; 
        % Moon when it was closer to the Earth 
        r_E_M=20; 
        RE=6370e3;
        Moon_cst.a=r_E_M*RE; 
        Moon_cst.e      = 0.0554;
        Moon_cst.i      = 5.16;  
        Moon_cst.raan   = 125.08;  
        Moon_cst.peri   = 318.15; 
        
        Moon_cst.obli=6.68;
        
        Moon_cst.hlid=50e3; 
        Moon_cst.hocean=50e3;
        Moon_cst.E=170e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=2800; 
        Moon_cst.name='Moon';
        Moon_cst.name_ab='Moon';

    % Martian moons
    case 41
        % Phobos
    case 42
        % Deimos
        
    % Jovian moons
    case 51
        % Io
        Moon_cst.Mp     = 1.8986E27; 
        Moon_cst.r      = 1821.46E3; 
        
        Moon_cst.obli=0.0021;
        Moon_cst.Mass=  8.9297e+22;
        Moon_cst.a      = 421800000;  
        Moon_cst.e      = 0.0041;
        Moon_cst.i      = 0.036;  
        Moon_cst.raan   = 43.977;  
        Moon_cst.peri   = 84.129;
        
        Moon_cst.hlid=30e3; 
        Moon_cst.hocean=50e3;
        Moon_cst.E=170e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=2800; 
        Moon_cst.name='Io';
        Moon_cst.name_ab='Io';
        
    case 52
        % Europa
        Moon_cst.Mp     = 1.8986E27; 
        Moon_cst.r      = 1560.8E3;
        
        
        Moon_cst.obli=0.053;
        Moon_cst.a      = 671100000;  
        Moon_cst.e      = 0.0094;
        Moon_cst.i      = 0.466;  
        Moon_cst.raan   = 219.106;  
        Moon_cst.peri   = 88.970;  
        
        
        Moon_cst.Mass   = 4.7998E22; 
        Moon_cst.C22    = 131.5E-6;
        Moon_cst.J2     = 10/3*Moon_cst.C22;
        Moon_cst.Im     = 0.3477;
        
        Moon_cst.hlid=(5+30)/2*1e3; 
        Moon_cst.hlid=30*1e3; 
        Moon_cst.hocean=(100+130)/2*1e3;
        Moon_cst.hocean=(100+100)/2*1e3;
        Moon_cst.E=8.8e9; 
        %Moon_cst.nu=0.33;
        Moon_cst.nu=0.499; % Europa proposal 
        Moon_cst.E=2*2e9*(1+Moon_cst.nu); % Europa proposal 
        Moon_cst.hlid=15*1e3; % Europa proposal 
        Moon_cst.hocean=123*1e3; % Europa proposal 
        Moon_cst.rho_o=1000; 
        Moon_cst.name='Europa';
        Moon_cst.name_ab='Eu';

    case 53
        % Ganymede
        Moon_cst.Mp     = 1.8986E27; 
        Moon_cst.r      = 2631.2E3;
        Moon_cst.Mass   = 14.8E22; 
        Moon_cst.a      = 1070400000;
        Moon_cst.e      = 0.0013;
        Moon_cst.i      = 0.177;
        Moon_cst.raan   = 63.552;
        Moon_cst.peri   = 192.417;
        Moon_cst.obli=0.0033;
        
        %Moon_cst.hlid=(20+120)/2*1e3; 
        Moon_cst.hlid=(125)*1e3; 
        Moon_cst.hocean=(150+150)/2*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000; 
        Moon_cst.name='Ganymede';
        Moon_cst.name_ab='Ga';
                    
    case 54
        % Callisto
        Moon_cst.Mp     = 1.8986E27; 
        Moon_cst.r      = 2410.E3;
        
        Moon_cst.obli=0.24;
        Moon_cst.Mass=10.8e22;
        Moon_cst.a      = 1882700000;
        Moon_cst.e      = 0.0074;
        Moon_cst.i      = 0.192;
        Moon_cst.raan   = 298.848;
        Moon_cst.peri   = 52.643;
        
        Moon_cst.hlid=(100+130)/2*1e3; 
        Moon_cst.hocean=(20+140)/2*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        Moon_cst.name='Callisto';
        Moon_cst.name_ab='Ca';
               
    % Saturnian moons
    case 61
        % Mimas
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 198.2E3;
        Moon_cst.Mass=3.7493e19;
        Moon_cst.a      = 185539000;
        Moon_cst.e      = 0.0196;
        Moon_cst.i      = 1.574;
        Moon_cst.obli=0.041;
        Moon_cst.raan   = 173.026;
        Moon_cst.peri   = 332.478;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        
        Moon_cst.hlid=(24+31)/2*1e3; 
        Moon_cst.hocean=71.6e3-Moon_cst.hlid;
        
        Moon_cst.name='Mimas';
        Moon_cst.name_ab='Mi';
        
        Moon_cst.mum=17e9/(2*(1+0.33));
        
               
    case 62
        % Enceladus
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 252.10E3;        
        Moon_cst.obli=0.00014;
        Moon_cst.obli=4.5e-4; % Matsuyama 2018
        Moon_cst.Mass=1.0799e+20;
        Moon_cst.a      = 238037000;
        Moon_cst.e      = 0.0047;
        Moon_cst.i      = 0.003;
        Moon_cst.raan   = 343.266;
        Moon_cst.peri   = 188.319;
        
        Moon_cst.hlid=(23+23)/2*1e3; 
        Moon_cst.hocean=(38+38)/2*1e3;
        Moon_cst.E=8.8e9; 
        %Moon_cst.E=9.31e9;
        Moon_cst.nu=0.33;
        Moon_cst.nu=0.49; %Europa proposal
        Moon_cst.E=2*2e9*(1+Moon_cst.nu); % Europa proposal 
        %Moon_cst.r      = 250E3; % Europa proposal 
        Moon_cst.rho_o=1000;
        Moon_cst.name='Enceladus';
        Moon_cst.name_ab='En';
        
        Moon_cst.mum=17e9/(2*(1+0.33));
                   
    case 63
        % Tethys
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 531.10E3;
        
        Moon_cst.a      = 294672000;
        Moon_cst.e      = 0.0001;
        Moon_cst.i      = 1.091;
        Moon_cst.raan   = 259.845;
        Moon_cst.peri   = 44.843;
        Moon_cst.name='Tethys';
        
        Moon_cst.density=0.98*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.obli=0.039;
        Moon_cst.hlid=(220)*1e3; 
        Moon_cst.hocean=(50)*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
                
    case 64
        % Dione
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 561.7E3;
        Moon_cst.obli=0.002;
        Moon_cst.Mass=1.0955e+21;
        Moon_cst.a      = 377415000;
        Moon_cst.e      = 0.0022;
        Moon_cst.i      = 0.028;
        Moon_cst.raan   = 290.615;
        Moon_cst.peri   = 284.111;
        
        Moon_cst.hlid=(99+99)/2*1e3; 
        Moon_cst.hocean=(65+65)/2*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        Moon_cst.name='Dione';
        Moon_cst.name_ab='Di';
        
        Moon_cst.mum=17e9/(2*(1+0.33));
               
    case 65
        % Rhea
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 764.4E3;
        
        Moon_cst.a      = 527068000;
        Moon_cst.e      = 0.001;
        Moon_cst.i      = 0.333;
        Moon_cst.raan   = 351.018;
        Moon_cst.peri   = 213.663;
        
        Moon_cst.density=1.24*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(401)*1e3; 
        Moon_cst.hocean=(16)*1e3;
        Moon_cst.obli=0.03;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        
                
    case 66
        % Titan
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 2574.73E3;
        
        Moon_cst.obli=0.32;
        Moon_cst.Mass=  1.3452e+23;
        Moon_cst.a      = 1221865000;
        Moon_cst.e      = 0.0288;
        Moon_cst.i      = 0.306;
        Moon_cst.raan   = 28.758;
        Moon_cst.peri   = 179.920;
        
        Moon_cst.hlid=(50+150)/2*1e3; 
        Moon_cst.hocean=(90+420)/2*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        Moon_cst.name='Titan';
        Moon_cst.name_ab='Ti';
                       
    case 67
        % Hyperion
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 133E3;
        
        Moon_cst.a      = 1500934000;
        Moon_cst.e      = 0.0232;
        Moon_cst.i      = 0.615;
        Moon_cst.raan   = 263.882;
        Moon_cst.peri   = 303.162;
              
    case 68
        % Iapetus
        Moon_cst.Mp     = 568.46E24;
        Moon_cst.r      = 735.6E3;
        
        Moon_cst.a      = 3560851000;
        Moon_cst.e      = 0.0293;
        Moon_cst.i      = 8.313;
        Moon_cst.raan   = 81.189;
        Moon_cst.peri   = 271.599;
               
    case 69
        % Phoebe 
        
    % Uranus' moons
    case 71
        % Miranda
        Moon_cst.Mp     = 86.832E24;
        Moon_cst.r      = 235.8E3;
        
        Moon_cst.a      = 129900000;
        Moon_cst.e      = 0.0014;
        Moon_cst.i      = 4.338;
        Moon_cst.raan   = 326.438;
        Moon_cst.peri   = 68.312;
        
        Moon_cst.density=1.18*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(115)*1e3; 
        Moon_cst.hocean=(10)*1e3;
        Moon_cst.obli=0.021;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;

    case 72
        % Ariel
        Moon_cst.Mp     = 86.832E24;
        Moon_cst.r      = 578.9E3;
        
        Moon_cst.a      = 190900000;
        Moon_cst.e      = 0.0012;
        Moon_cst.i      = 0.041;
        Moon_cst.raan   = 22.394;
        Moon_cst.peri   = 115.349;
        
        Moon_cst.density=1.54*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(190)*1e3; 
        Moon_cst.hocean=(19)*1e3;
        Moon_cst.obli=0.0005;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
               
    case 73
        % Umbriel
        Moon_cst.Mp     = 86.832E24;
        Moon_cst.r      = 584.7E3;
        
        Moon_cst.a      = 266000e3;
        Moon_cst.e      = 0.0039;
        Moon_cst.i      = 0.128;
        Moon_cst.raan   = 33.485;
        Moon_cst.peri   = 84.709;
        
        Moon_cst.density=1.52*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(190)*1e3; 
        Moon_cst.hocean=(10)*1e3;
        Moon_cst.obli=0.0026;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        
               
    case 74
        % Titania
        Moon_cst.Mp     = 86.832E24;
        Moon_cst.r      = 788.9E3;
        
        Moon_cst.a      = 436300000;
        Moon_cst.e      = 0.0012;
        Moon_cst.i      = 0.079;
        Moon_cst.raan   = 99.771;
        Moon_cst.peri   = 284.4;
        
        Moon_cst.density=1.65*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(230)*1e3; 
        Moon_cst.hocean=(20)*1e3;
        Moon_cst.obli=0.014;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        
               
    case 75
        % Obreron
        Moon_cst.Mp     = 86.832E24;
        Moon_cst.r      = 761.4E3;
        
        Moon_cst.a      = 583500000;
        Moon_cst.e      = 0.0014;
        Moon_cst.i      = 0.068;
        Moon_cst.raan   = 279.771;
        Moon_cst.peri   = 104.4;
        
        Moon_cst.density=1.66*1000;
        Moon_cst.Mass=4/3*pi*Moon_cst.density*Moon_cst.r^3;
        Moon_cst.hlid=(230)*1e3; 
        Moon_cst.hocean=(20)*1e3;
        Moon_cst.obli=0.075;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        
    % Neptune's moons
    case 81
        % Triton
        Moon_cst.Mp     = 102.43E24;
        Moon_cst.r      = 1353.E3;
        
        Moon_cst.obli=0.35;
        Moon_cst.Mass=  2.1390e+22;
        Moon_cst.a      = 354759000;
        Moon_cst.e      = 0.000016;
        Moon_cst.i      = 156.865;
        Moon_cst.raan   = 177.608;
        Moon_cst.peri   = 66.142;        
        
        Moon_cst.hlid=(150+150)/2*1e3; 
        Moon_cst.hocean=(150+150)/2*1e3;
        Moon_cst.E=8.8e9; 
        Moon_cst.nu=0.33;
        Moon_cst.rho_o=1000;
        Moon_cst.name='Triton';
        Moon_cst.name_ab='Tri';
    % Pluto's moons
    case 91
        % Charon
        Moon_cst.Mp     = 1.32E22;
        Moon_cst.r      = 603.6E3;
        
        Moon_cst.Mass=1.5327e+21;
        Moon_cst.a      = 17536000;
        Moon_cst.e      = 0.0022;
        Moon_cst.i      = 0.001;
        Moon_cst.raan   = 85.187;
        Moon_cst.peri   = 71.255;       
    case 92
        % Pluto
        Moon_cst.Mp     = 1.32E22;
        Moon_cst.r      = 1183.3E3;
        
        Moon_cst.Mass=1.3e+22;
        Moon_cst.a      = 17536000;
        Moon_cst.e      = 0.0022;
        Moon_cst.i      = 0.001;
        Moon_cst.raan   = 85.187;
        Moon_cst.peri   = 71.255;
        
    otherwise
        disp('selected moon is not included in the database');
end
Moon_cst.grav =Moon_cst.G*Moon_cst.Mass/Moon_cst.r^2;
Moon_cst.n      = sqrt(Moon_cst.G*Moon_cst.Mp/(Moon_cst.a^3));   
Moon_cst.Period = 2*pi/Moon_cst.n; 
Moon_cst.Ae=3.1963*(Moon_cst.n*Moon_cst.r)^2*Moon_cst.e/Moon_cst.grav;
Moon_cst.Aobli=1.4616*(Moon_cst.n*Moon_cst.r)^2*sind(Moon_cst.obli)/Moon_cst.grav;
Moon_cst.density=Moon_cst.Mass/(4/3*pi*Moon_cst.r^3);
Moon_cst.Lamb=Moon_cst.n^2*Moon_cst.r^2/(Moon_cst.grav*Moon_cst.hocean);
Moon_cst.rigidity=Moon_cst.E*Moon_cst.hlid/(Moon_cst.grav*Moon_cst.r^2*Moon_cst.rho_o); 
Moon_cst.ben_rigidity=1/(12*(1-Moon_cst.nu^2))*(Moon_cst.hlid/Moon_cst.r)^2;
%core radius
Moon_cst.rc=Moon_cst.r-Moon_cst.hlid-Moon_cst.hocean;
%core density 
Moon_cst.rho_c=(Moon_cst.r/Moon_cst.rc)^3*(Moon_cst.density-Moon_cst.rho_o*(1-(Moon_cst.rc/Moon_cst.r)^3));
%core gravity 
Moon_cst.gravc=Moon_cst.G*4/3*pi*Moon_cst.rho_c*Moon_cst.rc;
Moon_cst.muc=Moon_cst.mum/(Moon_cst.gravc*Moon_cst.rho_c*Moon_cst.rc);
end
