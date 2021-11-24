#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  excercise1.py
  

# EXPORTAR ANTES A PYTHON LA RUTA: PYTHONPATH=/home/esteban/CosmoMC_test/python:$PYTHONPATH
import getdist.plots as plots
g = plots.getSubplotPlotter(chain_dir = './chains/chains')#-ee-run-boost')


g.rectangle_plot(['omegabh2', 'omegach2', 'theta' , 'tau' , 'ns' , 'logA'  ], ['H0','omegam','sigma8'], 
   roots = [ 'boost_0_test_ee', 'boost_1_test_ee','boost_2_test_ee', 'boost_3_test_ee'], 
   filled=True, 
   legend_labels=['No boost', 'Nominal value','10 times Nominal value ', '100 times Nominal value'],
   legend_loc='upper right',
     );

g.export('./cmb_ee_only_boost_rectangle.png')


g.triangle_plot( ['boost_0_test_ee', 'boost_1_test_ee','boost_2_test_ee', 'boost_3_test_ee'], ['omegabh2', 'omegach2', 'theta' , 'ns' , 'H0' , 'omegamh2'],  
  filled=True, 
  legend_labels=['No boost', 'Nominal value','10 times Nominal value ', '100 times Nominal value'],
  legend_loc='upper right',
     );

g.export('./cmb_ee_only_boost_triangle.png')


