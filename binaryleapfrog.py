# %%
# import the packages
import numpy as np
import matplotlib.pyplot as plt

# %%
# the constants

G = 4*(np.pi**2) # gravitational constant units: AU^3/yr^2
M1 = M2 = 0.5 # solar mass ... makes equations easier
R = 1 # AU for circular orbits

# %%
# initialising variables / initial conditions

# initial position = x0 and y0 (in AU)
# initial velocity = v_x and v_y (in AU/yr)

# initial position
x0 = 1
y0 = 0

# initial velocity
v_x = 0
v_y = 4

# %%
# Using Kepler's Third Law, the equation of a period of a circular orbit is

T = np.sqrt(((4*(np.pi**2))/G*M1)*R**3)

# the time step - small iteration for the loop to show the approximate motion at that time for that x amount of time

dt = 0.0015
step = int(T/dt)    # how many years iteration

# to store the trajectories, use arrays
xvalues = [x0]
yvalues = [y0]

# initialise the variables so that it uses after the ones stored in the array
x = x0
y = y0

# to be used when the volecity is getting updated
vx = v_x
vy = v_y

energies = []
times = []

# %%
# position of the stars
star1 = np.array([-0.2, 0.0])   # left origin
star2 = np.array([0.2, 0.0])    # right origin

# Using the Leapfrog loop  

for i in range(step):
    # distance to star 1
    dx1 = x - star1[0]
    dy1 = y - star1[1]
    r1 = np.sqrt(dx1**2 + dy1**2)    # distance from the star circular orbit and modulus
    ax1 = -(G*M1*dx1)/(r1**3)    # acceleration in the x direction
    ay1 = -(G*M1*dy1)/(r1**3)    # acceleration in the y direction

# distance to star 2
    dx2 = x - star2[0]
    dy2 = y - star2[1]
    r2 = np.sqrt(dx2**2 + dy2**2)    # distance from the star circular orbit and modulus
    ax2 = -(G*M2*dx2)/(r2**3)          # acceleration in the x direction
    ay2 = -(G*M2*dy2)/(r2**3)          # acceleration in the y direction

    # total acceleration
    ax = ax1+ax2
    ay = ay1+ay2

    # half-step velocity
    vx_half = vx + 0.5*ax*dt
    vy_half = vy + 0.5*ay*dt

    # update position
    x = x + vx_half*dt
    y = y + vy_half*dt

    # updated acceleration
    dx1 = x - star1[0]
    dy1 = y - star1[1]
    r1 = np.sqrt(dx1**2 + dy1**2)
    ax1 = -(G*M1*dx1)/(r1**3)    # acceleration in the x direction
    ay1 = -(G*M1*dy1)/(r1**3)    # acceleration in the y direction


    dx2 = x - star2[0]
    dy2 = y - star2[1]
    r2 = np.sqrt(dx2**2 + dy2**2)
    ax2 = -(G*M2*dx2)/(r2**3)          # acceleration in the x direction
    ay2 = -(G*M2*dy2)/(r2**3)          # acceleration in the y direction

    # total updated acceleration
    ax_new = ax1+ax2
    ay_new = ay1+ay2

    # full half-step new velocity
    vx = vx_half + 0.5*ax_new*dt
    vy = vy_half + 0.5*ay_new*dt

    # adds the new trajectories to the end of the list
    # snapshot interval to only add every 2nd value
    snap = 20
    if i % snap == 0:
        xvalues.append(x)
        yvalues.append(y)

    # energy time graph
    KE = 0.5*(M1 + M2)*(vx**2 + vy**2)
    PE = -G*(M1)/r1 - G*(M2)/r2

    energies.append(KE + PE)
    times.append(i*dt)
    

# %%
plt.figure(figsize=(5.5,14))

# plotting the loop
plt.subplot(3,1,1)
plt.plot(xvalues, yvalues, 'x', markersize = 3, label = "Celestial Object") # orbit of the object
plt.plot(star1[0], star1[1], "o", markersize = 10, label = "Star 1")
plt.plot(star2[0], star2[1], "ro", markersize = 10, label = "Star 2")
plt.title('Binary Orbital System with Leapfrog Method')
plt.xlabel("x (au)")
plt.ylabel("y (au)")
plt.legend(loc='upper left')
plt.grid()

# energy time graph
plt.subplot(3,1,2)
plt.plot(times, energies, color='red')
plt.ylim(min(energies)-0.1, max(energies)+0.1)
plt.title('Orbit Energy vs. Time with RK2')
plt.xlabel("Time (yr)")
plt.ylabel("Total Energy (AU$^2$/yr$^2$)")
plt.grid()

# energy time plot zoomed in
plt.subplot(3,1,3)
plt.plot(times, energies, color='red')
plt.title('Circular Orbit Energy vs. Time with Leapfrog Zoomed In')
plt.xlabel("Time (yr)")
plt.ylabel("Total Energy (AU$^2$/yr$^2$)")
plt.grid()

plt.tight_layout()

# %%
plt.figure(figsize=(14,6))
plt.rcParams.update({'font.size': 11})

# plotting the loop
plt.subplot(1,2,1)
plt.plot(xvalues, yvalues, 'x', markersize = 3, label = "Celestial Object") # orbit of the object
plt.plot(star1[0], star1[1], "o", markersize = 10, label = "Star 1")
plt.plot(star2[0], star2[1], "ro", markersize = 10, label = "Star 2")
plt.title('Binary Elliptical Orbital System (Leapfrog Method)')
plt.xlabel("x (au)")
plt.ylabel("y (au)")
plt.legend(loc='upper right')
plt.grid()

# energy time graph
plt.subplot(1,2,2)
plt.plot(times, energies, color='red')
plt.ylim(min(energies)-0.1, max(energies)+0.1)
plt.title('Binary Elliptical Orbital System Energy vs. Time (Leapfrog Method)')
plt.xlabel("Time (yr)")
plt.ylabel("Total Energy (AU$^2$/yr$^2$)")
plt.grid()

plt.tight_layout()


