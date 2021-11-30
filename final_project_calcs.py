import math
import numpy as np
import pylab as pl

print("\n------ SIMULATION INITIATING ------\n")
pl.close('all') 

# Planetary Constants
mu = 5.788 * 10**6 # gravitational parameter, km^3/sec^2
SG_coef = 0.6645 * 10**-4 # Sutton Graves coefficient for Uranus, kg^.5/m
g_earth = 9.81 # earth gravity, m/s^2
R = 3164.92 # gas constant for atmospheric constituents, J/kg K
gamma = 1.4 # ratio of specific heats, assumed for H2 majority atmosphere
T_avg = 78.15 # avg temperature, K
a_avg = math.sqrt(gamma * R * T_avg) # average speed of sound, m/s
rho_0 = 0.42 # density at 1 Bar level, kg/m^3
r_uranus = 24973 # equatorial radius to 1 Bar level, km
H_uranus = 2.77 * 10**4 # scale height, m
H_earth = 8500 # scale height, m
A_earth = 1/H_earth # 1/m
A_uranus = 1/H_uranus # 1/m
vk_p = 101325 * math.exp(-A_earth * 100000) # estimated pressure at Earth's Von Karman line, Pa
vk_h_uranus = math.log(vk_p /101325)/(-A_uranus)/1000 # Von Karman line altitude on Uranus, km

print("Pressure at Von Karman Line: " + str('%.3f' % vk_p) + " Pa")
print("Von Karman Line Altitude Equivalent for Uranus: " + str('%.3f' % vk_h_uranus) + " km above 1 Bar altitude\n")

# Mothership Constants & Velocity Calculations
h_mothership = 500 # altitude of mothership above karman line, km
r_a = r_uranus + h_mothership # radius of apoapsis, km
r_p = r_uranus # radius of periapsis, km
semi_maj_axis = (r_a + r_p)/2 # semi-major axis, km
v_mothership = math.sqrt(mu/(r_a)) # mothership orbital velocity (circular), km/s
v_a = math.sqrt((2*mu/r_a) - (mu/semi_maj_axis)) # apoapsis veloicty based on hohmann deorbit, km/s
v_p = math.sqrt((2*mu/r_p) - (mu/semi_maj_axis)) # periapsis velocity - this would be the highest velocity encountered, km/s
dv = v_mothership - v_a # delta_v required to deorbit from mothership

print("Mothership Orbital Velocity: " + str('%.3f' % v_mothership) + " km/s @ altitude of " + str(h_mothership) +" km")
print("Velocity @ Apoapsis for deorbit: " + str('%.3f' % v_a) + " km/s")
print("Max Velocity from Hohmann deorbit: " + str('%.3f' % v_p) + " km/s")
print("Delta V required for deorbit: " + str('%.3f' % dv) + " km/\n")

# Probe Constants
c_d = 1.8 # coefficient of drag
d_probe = 1.5 # probe diameter, m
r_n = d_probe/4 # probe nose radius, m
m_probe = 660 # probe mass, kg
S_probe = math.pi/4 * d_probe**2 # probe surface area, m^2
beta = m_probe/(c_d * S_probe) # ballistic coefficient, kg/m^2
print("Ballistic Coefficient: " + str('%.2f' % beta) + " kg/m^2")
print(S_probe)

# EDL Calculations
v_e = v_p # assume that entry velocity is veloicty at periapsis, km/s
FPA_deg = np.array((-20, -23, -30, -35)) # varying FPA for further analysis, deg
FPA_rad = FPA_deg / 180 * math.pi

v_q_max = 0.846 * v_e # velocity at peak heating, km/s
v_n_max = 0.606 * v_e # velocity at peak deceleration, km/s
print("Velocity at peak heating: " + str('%.2f' % v_q_max) + " km/s")
print("Velocity at peak deceleration: " + str('%.2f' % v_n_max) + " km/s\n")

# Intialize Vectors
alt_vec = np.linspace(h_mothership,0,10001) # altitude vector based on mothership height, km
v_curr = np.zeros((len(alt_vec), len(FPA_deg))) # velocity vector
a_curr = np.zeros((len(alt_vec), len(FPA_deg))) # acceleration vector
a_curr_g = np.zeros((len(alt_vec), len(FPA_deg))) # acceleration vector, in g's
M_curr = np.zeros((len(alt_vec), len(FPA_deg))) # Mach vector
q_flux = np.zeros((len(alt_vec), len(FPA_deg))) # heat flux vector
del_s = np.zeros((len(alt_vec), len(FPA_deg))) # downrange distance vector
press = np.zeros((len(alt_vec), len(FPA_deg))) # dynamic pressure vector

for j in range(len(FPA_deg)):
    print("---------- FPA = " + str('%.1f' % FPA_deg[j]) + " deg ----------")
    Q_total = (SG_coef * (v_e*1000)**2 * math.sqrt((-beta * math.pi)/(r_n * A_uranus * math.sin(FPA_rad[j]))))/100**2 # total heat load, J/cm^2
    m_TPS = 0.091 * Q_total**0.51575 # TPS mass fraction based on Lecture 5 notes, %
    q_max = (SG_coef * (0.6055 * ((v_e*1000)**3)) * math.sqrt((-1/r_n) *(beta * A_uranus * math.sin(FPA_rad[j]) / 3)))/1000 # max heat flux, kW/m^2
    q_max_cm = q_max * 1000 / (100**2) # convert from kW/m^2 to W/cm^2
    n_max = ((v_e*1000)**2 * A_uranus * math.sin(-FPA_rad[j]))/(2 * g_earth * math.exp(1)) # peak deceleration, g's
    h_n_max = (1 / A_uranus * math.log(-rho_0 / (beta * A_uranus * math.sin(FPA_rad[j]))))/1000 # altitude of peak deceleration, km
    
    C = rho_0/(2 * beta * A_uranus * math.sin(FPA_rad[j])) # arbitrary C coefficient for velocity calc
    h_q_max = (math.log(math.log((v_q_max) / (v_e)) / C) / (-A_uranus))/1000 # altitude of peak heating, km
    
    print("Total heat load: " + str('%.2f' % Q_total) + " J/cm^2")
    print("TPS mass fraction percentage: " + str('%.2f' % m_TPS) + "%")
    print("Peak heating: " + str('%.2f' % q_max) + " kW/m^2")
    print("Peak heating in W/cm^2: " + str('%.2f' % q_max_cm) + " W/cm^2")
    print("Altitude of peak heating: " + str('%.2f' % h_q_max) + " km")
    print("Peak deceleration: " + str('%.2f' % n_max) + " g's")
    print("Altitude of peak deceleration: " + str('%.2f' % h_n_max) + " km")
    #print("Ballistic Downrange Distance: " + str('%.2f' % del_s) + " km")
    
    flag = 0
    for i in range(int(len(alt_vec))):
        v_curr[i,j] = ((v_e*1000) * math.exp(C * math.exp(-A_uranus * alt_vec[i]*1000)))/1000 # velocity calculation for ballistic entry, km/s
        M_curr[i,j] = v_curr[i,j] * 1000 / a_avg # velocity in terms of Mach
        a_curr[i,j] = (-v_e*1000)**2 * math.sin(FPA_rad[j]) * C * A_uranus * math.exp(2 * C * math.exp(-A_uranus * alt_vec[i]*1000)) * math.exp(-A_uranus * alt_vec[i]*1000) # deceleration, m/s^2
        a_curr_g[i,j] = a_curr[i,j] / g_earth # deceleration, g's                                                                      
        q_flux[i,j] = (SG_coef * ((v_curr[i,j]*1000)**3) * math.sqrt((rho_0 * math.exp(-A_uranus * alt_vec[i]*1000)) / r_n) )/(100**2) # heat flux as a factor of altitude & velocity, W/cm^2
        del_s[i,j] = (math.cos(-FPA_rad[j]) / math.sin(-FPA_rad[j])) * (400 - alt_vec[i]) # downrange distance over time, km
        press[i,j] = (0.5 * rho_0 * math.exp(-A_uranus * alt_vec[i]*1000) * (v_curr[i,j]*1000)**2)/1000 # dynamic pressure as a functional of descent, kPa 
        
        if flag == 0:
            if M_curr[i,j] < 0.9:
                print("Altitude when Mach < 0.9: " + str('%.1f' % alt_vec[i]) + " km\n")
                #print("Velocity when Mach < 0.9: " + str('%.2f' % (v_curr[i,j]*1000)) + " m/s\n")
                flag = 1
    
    #title_1 = 'Altitude vs Velocity - FPA=' + str(FPA_deg[j]) + ' deg'
    #title_2 = 'Altitude vs Heat Flux - FPA=' + str(FPA_deg[j]) + ' deg'
    
# generate plots
pl.figure()       
pl.plot(v_curr, alt_vec)
pl.xlabel('Velocity (km/s)')
pl.ylabel('Height (km)')
pl.title('Altitude vs Velocity as a function of FPA')
pl.legend(FPA_deg)

pl.figure()       
pl.plot(a_curr_g, alt_vec)
pl.xlabel('Deceleration (g\'s)')
pl.ylabel('Height (km)')
pl.title('Altitude vs Deceleration as a function of FPA')
pl.legend(FPA_deg)

# pl.figure()       
# pl.plot(press, alt_vec)
# pl.xlabel('Dynamic Pressure (kPa)')
# pl.ylabel('Height (km)')
# pl.title('Altitude vs Dynamic Pressure as a function of FPA')
# pl.legend(FPA_deg)

# pl.figure()
# pl.plot(M_curr, alt_vec)
# pl.xlabel('Mach Number')
# pl.ylabel('Height (km)')
# pl.title('Altitude vs Mach Number as a function of FPA')
# pl.legend(FPA_deg)

pl.figure()
pl.plot(q_flux, alt_vec)
pl.xlabel('Heat Flux (W/cm^2)')
pl.ylabel('Height (km)')
pl.title('Altitude vs Heat Flux as a function of FPA')
pl.legend(FPA_deg)

pl.figure()
pl.plot(del_s, alt_vec)
pl.xlabel('Downrange Distance (km)')
pl.ylabel('Height (km)')
pl.title('Altitude vs Downrange Distance as a function of FPA')
pl.legend(FPA_deg)
pl.xlim(0)
pl.ylim(0)

print("------ MISSION SUCCESS! FLOATING THROUGH THE SKIES OF URANUS! ------")



