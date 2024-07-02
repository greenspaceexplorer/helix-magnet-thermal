# HELIX magnet thermal model class
# Author: Noah Green
# Last Modified: 7/18/2016

import sys
import math
import matplotlib
import matplotlib.pyplot as plt
from configobj import ConfigObj
#
#class ThermalSystem(object):
#    def __init__(self,fileout,file_materials,file_elements,file_interface):
#        
#        
#        
#
#        
#
#class Interface(object):
#    def __init__(self,element1,element2,intfclist):
#        """
#        Instantiator for thermal interface
#        Input: element, element, list(int)
#        note: for type, conduction = 1, convection = 2, radiation = 3
#        """
#        

###############################################################################
###############################################################################

class Material(object):
    def __init__(self,matlist):
        """
        Instantiator for type of material
        Input: list(string, number   , number    , number  , number , int  , number      )
        i.e.   list(name  , thermcond, emissivity, specheat, density, state, vaporization)
        note: for state, solid = 1, liquid = 2, gas = 3
        """
        
        self.__lst__ = matlist
        self.__name__ = matlist[0]
        self.__thermcond__ = float(matlist[1])
        self.__emissivity__ = float(matlist[2])
        self.__specheat__ = float(matlist[3])
        self.__density__ = float(matlist[4])
        self.__state__ = int(matlist[5])
        self.__vap__ = float(matlist[6])
    
        self.__valid__ = False
        self.IsValid()
        
        if not self.__valid__:
            raise ValueError('Invalid material list values')

        
    def __str__(self):
        name = "Name: {}".format(self.__name__)
        thermcond = "Thermal Conductivity: {}".format(self.__thermcond__)
        emissivity = "Emissivity: {}".format(self.__emissivity__)
        specheat = "Specific Heat: {}".format(self.__specheat__)
        density = "Density: {}".format(self.__density__)
        state = "State: {}".format(self.__state__)
        vap = "Heat of Vaporization: {}".format(self.__vap__)
        
        return "{}\n{}\n{}\n{}\n{}\n{}\n{}".format(name,thermcond,\
        emissivity,specheat,density,state,vap)

       
    def IsValid(self):
        """
        Checks to see if material property list is valid
        """
        self.__valid__ = True
        
        if type(self.__name__) != str:
            self.__valid__ = False
            print("Invalid name")
    # Setting functions
    def SetThermCond(self,tc):
        self.__thermcond__ = float(tc)
    def SetEmissivity(self,emis):
        self.__emissivity__ = float(emis)
    def SetDensity(self,dens):
        self.__density__ = float(dens)

    # Getting functions        
    def Name(self):
        return self.__name__
    def ThermCond(self):
        return self.__thermcond__
    def Emissivity(self):
        return self.__emissivity__
    def SpecHeat(self):
        return self.__specheat__
    def Density(self):
        return self.__density__
    def State(self):
        return self.__state__
    def Vaporization(self):
        return self.__vap__
    def PropertyList(self):
        return self.__lst__

            
class Element(Material):
    def __init__(self,mat,temp,mass,volume):
        """
        Instantiator for thermal element
        Input: Material, float, float, float, float 
        """
        
        Material.__init__(mat.property_list())
        self.__T__ = float(temp)
        self.__M__ = float(mass)
        self.__V__ = float(volume)
        if self.__M__ == -1 and self.__V__ != -1:
            self.__M__ = self.__V__*self.density()
        self.__ThCap__ = self.__M__*self.specheat()
    def __str__(self):
        mat_str = "Material: " + self.name()
        temp_str = "Temp(K): " + str(self.__T__)
        mass_str = "Mass(kg): " + str(self.__M__)
        vol_str = "Volume(m3): " + str(self.__V__)
        return "{}\n{}\n{}\n{}\n".format(mat_str,temp_str,mass_str,vol_str)
    # Setting functions
    def SetT(self,temp):
        self.__T__ = float(temp)
    def SetV(self,vol):
        self.__V__ = float(vol)
 #   def SetM
        