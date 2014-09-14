"""configuration file"""
import os
import pdb
import scipy as SP

CFG = {}

colors = {}
colors['cc'] = 'Peru'
colors['bio'] = 'DarkMagenta'
colors['tech'] = '#92c5de'
colors['Th2'] = 'DarkBlue'
colors['Th2_cc'] = 'RoyalBlue'
legends = {}
legends['cc'] = 'cell cycle'
legends['bio'] = 'biol. var'
legends['tech'] = 'tech. var'
legends['Th2'] = 'Th2 diff.'
legends['Th2_cc'] = 'interaction'

#color definitions and legends
CFG['colors'] = colors
CFG['legends'] = legends
#definitions for variance components
CFG['var_comp_fields_order'] = ['cc','bio','tech','Th2','Th2_cc']
CFG['var_comp_fields']  = SP.array([[CFG['var_comp_fields_order'].index(d),legends[d],colors[d]] for d in CFG['var_comp_fields_order']],dtype='object')
