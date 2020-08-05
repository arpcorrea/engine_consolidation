# -*- coding: utf-8 -*-
import yaml, sys
import os.path

class Deck():

    def __init__(self, inputhpath):
        if not os.path.exists(inputhpath):
            print("File " + inputhpath)
            sys.exit(1)
        else:
            with open(inputhpath,'r') as f:
                ## Container of the tags parsed from the yaml file
                self.doc = yaml.load(f, Loader=yaml.BaseLoader)
        self.create_folder_structure()
        self.type=self.doc["Problem Type"]["Name"]
        
        
    def create_folder_structure(self):        
        self.plot_dirTemp = self.doc["Plot"]["Temp Output Folder"]
        self.plot_dirDic = self.doc["Plot"]["Dic Output Folder"]
        check_folder_Temp = os.path.isdir(self.plot_dirTemp)
        check_folder_Dic = os.path.isdir(self.plot_dirDic)
        if not check_folder_Temp:
              os.makedirs(self.plot_dirTemp)
        if not check_folder_Dic:       
              os.makedirs(self.plot_dirDic)

