# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:16:38 2016

@author: green
"""
from configobj import ConfigObj

def main():
    filename = "Materials.dat"
    config = ConfigObj(filename)
    for key in config:
        print(key)
        subdict=config[key]
        for subkey in subdict:
            print(subkey)
            
    
main()