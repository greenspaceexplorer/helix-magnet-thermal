from ThermalClass import *


def main():
    """
    Main function for HEAT magnet thermal calculation
    """
    name = 'aluminum'
    thermcond = 117 # W/m-K
    emissivity = 0.06
    specheat = 900 # J/kg-K
    density = 2660 # kg/m3
    state = 1
    vap = -1
    
     
    al_list = [name,thermcond,emissivity,specheat,density,state,vap]
 
    aluminum = Material(al_list)
     
    print(aluminum)



main()
