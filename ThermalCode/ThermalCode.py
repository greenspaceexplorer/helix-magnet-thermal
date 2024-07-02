# Thermal model code for HELIX magnet
# Author: Noah Green
# Translated from QBasic code in Patrick Koehn's thesis

#from __future__ import print_function
import math
import sys
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt

def make_plot(x_list, y_list, title_plot, title_x, title_y, save_name):
    """
    Makes plot with given coordinates
    Input: list, list, string, string, string, string
    Output: void
    """
    
    # Get plot boundaries
    x_min = min(x_list)
    x_max = max(x_list)
    y_min = min(y_list) - 0.1 * abs(min( y_list ))
    y_max = max(y_list) + 0.1 * abs(max( y_list ))
    
    # Plot points
    figure, (ax1) = plt.subplots(1,1)
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    ax1.set_xlabel(title_x)
    ax1.set_ylabel(title_y)
    ax1.set_title(title_plot)
    ax1.plot(x_list, y_list)
    figure.savefig(save_name)
    plt.close(figure)

def height(x):
    """
    Returns height of LHe for a given volume
    input: float
    output: float
    """
    if x <= 43 and x >= 0:
        return (0.231166897578547 + 0.19579947200452*x\
            - 0.00147163230560322*x**2)
        
    elif x > 43 and x < 231:
        return ( -22.4033436605957 + 1.25361724755983*x\
            -  0.0183301586082401*x**2 + 0.000138923728362529*x**3\
            -  5.16376654974849e-7*x**4 + 7.53834532810923e-10*x**5)
            
    elif x >= 231 and x <= 277:
        return (104.127497960411 - 0.708174149227273*x\
            + 0.00163976643013219*x**2)
        
    else:
        sys.exit("Error: invalid volume")

# todo: normalize this to actual contact area
def contact_area(x):
    """
    Returns approximation of the area (in^2) of LHe in contact with wall for 
        given volume 
    input: float
    output: float
    """
    if x <= 43 and x >= 0:    
        #return integrate.quad(ca_low,1.e-13,volume)
        return 13.3183289108*x + 0.8507808080685*x**2 - 0.02866324423394*x**3\
            + 0.000435708787407*x**4 - 2.966784840519e-6*x**5
        
    elif x > 43 and x < 231:
        #return (integrate.quad(ca_low,1.e-13,43) \
        #    + integrate.quad(ca_mid,43,volume))
        return  -6979.4556826581 + 334.59508916278*x - 4.8875231024782*x**2\
            + 0.036912724297117*x**3 - 0.00013671140727553*x**4\
            + 1.9884834242872e-7*x**5
    
    elif x >= 231 and x <= 277:
        #return (integrate.quad(ca_low,1.e-13,43) \
        #    + integrate.quad(ca_mid,43,231)\
        #    + integrate.quad(ca_high,231,volume))
        return 5.4233262998781e-1*x**2 - 2.3216422554555e2*x\
            + 3.0714876200494e4
    
    else:
        sys.exit("Error: invalid volume")
        
    
    
def ca_low(x):
    """
    integrand for contact area function
    """
    a = math.sqrt(19)/9
    b = math.sqrt(12.2)/5.6
    c = math.sqrt(15.93)/6.35
    w1 = 10.5
    w2 = 2
    l = 60
    vmax = 277
    hmax = height(vmax)/2
    h = height(x)
    k1 = 2*math.sqrt(1+a**2) + math.sqrt(1+b**2) + math.sqrt(1+c**2)
    k2 = 2*a + b + c
    k3 = w1 + w2        
    
    dhdv = 0.19579947200452 - 0.00294326461120644*x
    return 2*(k1*math.sqrt(2*hmax*h - h**2) + k2*h + k3) * dhdv
 

def ca_mid(x):
    a = math.sqrt(19)/9
    b = math.sqrt(12.2)/5.6
    c = math.sqrt(15.93)/6.35
    w1 = 10.5
    w2 = 2
    l = 60
    vmax = 277
    hmax = height(vmax)/2
    h = height(x)
    k1 = 2*math.sqrt(1+a**2) + math.sqrt(1+b**2) + math.sqrt(1+c**2)
    k2 = 2*a + b + c
    k3 = w1 + w2  

    dhdv = 1.25361724755983 - 0.0366603172164802*x\
        -  0.000416771185087587*x**2 - 2.065506619899396e-6*x**3\
        -  3.769172664054615e-9*x**4
    return 2*(2*math.sqrt(2*hmax*h - h**2)\
        + l*hmax/math.sqrt(2*hmax*h - h**2))*dhdv
      
def ca_high(x):
    a = math.sqrt(19)/9
    b = math.sqrt(12.2)/5.6
    c = math.sqrt(15.93)/6.35
    w1 = 10.5
    w2 = 2
    l = 60
    vmax = 277
    hmax = height(vmax)/2
    h = height(x)
    k1 = 2*math.sqrt(1+a**2) + math.sqrt(1+b**2) + math.sqrt(1+c**2)
    k2 = 2*a + b + c
    k3 = w1 + w2        
    
    dhdv = 0.0032*x - 0.7082
    return 2*(k1*math.sqrt(2*hmax*h - h**2) + k2*(hmax-h) + k3) * dhdv

def warmup(Q1,THe,initperc,r,m,Helevel,perc,boil,expand,\
           marker,dV,flag,Qeach,block,LHeT,moles):
    """
    Warmup function for LHe < bp
    input: 11 int/floats, 5 lists
    output: 4-tuple
    """

    #level of LHe in the cryostat in inches
    #ht = 5.6834e-6 * Helevel**3 - 0.0023 * Helevel**2 + 0.358 * Helevel - 3.881
    #fraction = 1. - math.acos( ht/r - 1 )/math.pi
    fraction = 1. - math.acos( height(Helevel)/r - 1 )/math.pi
    #fraction = contact_area(Helevel)/contact_area(277)
#1. - math.atan(math.sqrt((r**2/(ht - r)**2) - 1.))/math.pi
    Qavail = Q1 * fraction
    perc = 0.
    dV = 0.

    for n in range(m):
        #skip if block is at bp
        if flag[n] == 1:
            continue
        #partition energy
        Qeach[n] = block[n] * Qavail / Helevel
        if n == (m-1):
            Qeach[n] += Q1 - Qavail
        #LHe heat capacity
        capacity = 3.1 * LHeT[n] + 1.5
        #change in temperature
        LdT = Qeach[n] / (moles[n] * capacity)
        LHeT[n] += LdT
        #new block size from thermal expansion
        block[n] *= (1. + 3. * expand * LdT)
        #increase in LHe volume
        dV += block[n] * 3. * expand * LdT
        #new percentage
        perc += block[n]/Helevel 
        if LHeT[n] >= THe:
            #flag if block is at bp
            flag[n] = 1
            #advance marker
            marker += 1
            #adjust perc to previous level
            perc -= block[n]/Helevel

    #adjust for roundoff error
    if marker == 0:
        perc = initperc
    #reset perc
    if marker == m:
        perc = 0;

    boil = (Qavail * (1. - perc) + (Q1 - Qavail) * flag[m-1]) * 0.000383
    return (perc,boil,marker,dV)
    
def main():
    """
    Main function of thermal model
    """ 
    ###########################################################################
    # User Inputs
    ###########################################################################

    #output file
    filename = "thermaloutput.csv"
    #initial temperature of shield
    T = 81.15
    #initial level of LHe in cryostat
    Helevel = 265
    #fraction of total boiloff passing through shield tube
    f = .5
    #shield heat-up error factor
    d = 0.45
    #initial percentage of LHe below bp
    initperc = 1.0
    #lower bound on temperature strats
    lower = 4.12
    #number of bins for temperature strats
    m = 10
    #ambient atmospheric pressure
    atmos = 14.504 #12.959
    #cryostat pressure
    stack = 3.5
    #room temperature
    Tskin = 293.
    #LHe volume where magnet quenches
    cutoff = 40.

    ###########################################################################
    # Fixed Constants
    ###########################################################################

    #radius of cryostat
    r = 17
    #initialize deltaT and time
    deltaT = 0.
    minutes = 0
    seconds = 0

    marker = 0
    if initperc == 0:
        marker = m

    ###########################################################################
    # Initial Calculations
    ###########################################################################

    #save initial Helevel
    Heinit = Helevel
    #initialize percent He below bp
    perc = initperc
    #calculate LHe bp from pressure
    THe = 0.069 * (atmos + stack) + 3.2
    #calculate expansion coefficient
    expand = 0.0624 * THe - 0.1243
    #calculate upper bound for temperature strats
    higher = THe - 0.02
    #calculate constant heat flux from support
    Qsupport = (Tskin - THe) * 0.001254
    #calculate constant heat flux from bore
    Qthroat = (Tskin**4 - THe**4) * 3.97e-11


    #Setup Temperature Proile Below Boiling Point
    mult = list(1./m for x in range(m))

    #Setup Temperature Distribution Below Boiling Point

    LHeT = list()
    flag = list()
    Qeach = list()
    block = list()
    moles = list()
    inc = 0

    #fill lists
    for n in range(m):
        LHeT.append(lower+inc)
        inc += (higher - lower)/m
        block.append(mult[n] * initperc * Heinit)
        moles.append((39.745 * LHeT[n] - 1.317) * block[n])
        flag.append(0)
        Qeach.append(0)

    fileout = open(filename,"w")
    print("{},{},{},{},{},{}".format("seconds","Helevel","Temp","boiloff",\
        "perc","marker"),file = fileout)

    ###########################################################################
    # Main Program
    ###########################################################################


    #initialize variables outside loop for proper scope
    Q1 = 0.
    boil = 0.
    level_list = []
    temp_list = []
    time_list = []
    marker_list = []
    boil_list = []
    perc_list = []
    while Helevel > cutoff:
        seconds += 1
        T += deltaT
        Qshield = 8.612e-9 * (T**4 - THe**4)
        Q1 = Qsupport + Qthroat + Qshield
        dV = 0.
        if marker < m:
            warmup_list = warmup(Q1,THe,initperc,r,m,Helevel,perc,boil,expand,\
                                 marker,dV,flag,Qeach,block,LHeT,moles)
            perc = warmup_list[0]
            boil = warmup_list[1]
            marker = warmup_list[2]
            dV = warmup_list[3]
        else:
            boil = Q1 * 0.000383
        
        Qinsul = 1292.06 * (1.e-12 * (Tskin**4 - T**4) + 4.2e-5 * (Tskin - T))
        coolingtube = (559. * T - 2082.81) * boil * f
        Alcap = -1.77e-7 * T**3 + 0.0000377 * T**2 + 0.004 * T - 0.129
        deltaT = d * (Qinsul - Qshield - coolingtube)/(19084.6075 * Alcap)

        Helevel += dV - boil
        
        if marker == m:
            dV = 0
        
        # Save values
        print("{},{},{},{},{},{}".format(seconds,Helevel,T,boil,perc,marker),\
          file = fileout)
        #if seconds%3600 == 0:
        #    print(seconds/3600, " hours")
        if seconds%1 == 0:
            time_list.append(seconds/3600)
            level_list.append(Helevel)
            temp_list.append(T)
            marker_list.append(marker)
            boil_list.append(boil*3600)
            perc_list.append(perc)

    ###########################################################################
    # Make output and close file(s)
    ###########################################################################
   
    make_plot(time_list, temp_list, "Shield Temperature", "Time (hrs)"\
            , "Temperature (K)", "shield_temp.png")
    make_plot(time_list, level_list, "LHe Volume", "Time (hrs)", "Volume (l)"\
            , "volume.png")
    make_plot(time_list, marker_list, "Marker Value", "Time(hrs)", "Marker"\
            , "marker.png")
    make_plot(time_list, boil_list, "LHe Boiloff Rate","Time(hrs)"\
            , "Rate(l/hr)","boil.png")
    fileout.close()
    print("Data saved in {}".format(filename))
    print("Hold time: {} days, {} hours, {} minutes\
    , {} seconds ( or {} hours )".format(seconds//86400\
            ,(seconds%86400)//3600,((seconds%86400)%3600)//60\
                ,((seconds%86400)%3600)%60,seconds/3600))
    

main()

