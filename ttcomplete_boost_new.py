#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# scales.py
#	This code aims to analyse chains from COSMOMC sampler using getdist libraries
#	Be sure that you have installed all the necessary libraries
#	The code produces statistical plots used to check convergence in Markov Chains
#	Gelman-Rubin test is calculated for the convergence baseline of 0.01  this test takes longer 
#	in calculating the convergence 1-R 
#
#
# EXPORTAR ANTES A PYTHON LA RUTA: PYTHONPATH=/home/esteban/CosmoMC_test/python:$PYTHONPATH

import math
import numpy as np
import getdist.plots as plots
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib
import scipy 
import pandas as pd
import subprocess
import sys
import scipy.stats

from scipy.stats import norm
from matplotlib import rc
from getdist import loadMCSamples
from getdist import loadMCSamples
from getdist import covmat
from getdist import MCSamples
from tabulate import tabulate
from scipy.optimize import curve_fit
from matplotlib.projections.geo import GeoAxes
from mpl_toolkits.mplot3d import Axes3D


def font_stile(a):	
	rc('font',**{'family':'serif','serif':[a]})
	rc('text', usetex=True)
	#matplotlib.rcParams.update({'font.size': 18})
	return True
	
def correlations_plot_fun(folder, root0, root1, root2, label0, label1,  label2, output_name, flag, true):
	if flag=='Triangle_Harmonic':
		g = plots.getSubplotPlotter(chain_dir = folder)
		g.settings.figure_legend_frame = True
		g.settings.alpha_filled_add = 0.4
		g.settings.title_limit_fontsize = 14
		g.triangle_plot([ root0, root1, root2],['beta10', 'beta11r', 'beta11i'], 
			markers={'beta10':0.0},
			filled=True, 
			legend_labels=[label0, label1, label2], 
			legend_loc='upper right',
			contour_colors=['red','green','darkblue'],
			title_limit=1)
		n = ['beta10', 'beta11r', 'beta11i']
		h0 = loadMCSamples('./moremcmc_new/'+root0)
		true1 = h0.mean(n, where=None)
	
		h1 = loadMCSamples('./moremcmc_new/'+root1)
		true3 = h1.mean(n, where=None)
		
		h2 = loadMCSamples('./moremcmc_new/'+root2)
		true2 = h2.mean(n, where=None)
	
	
		for ax in g.subplots[:,0]:
			ax.axvline(true[0], color='black', ls='--')
			ax.axvline(true1[0], color='red', ls='--')
			ax.axvline(true2[0], color='blue', ls='--')
			ax.axvline(true3[0], color='lime', ls='--')
		for ax in g.subplots[1:,1]:
			ax.axvline(true[1], color='black', ls='--')
			ax.axvline(true1[1], color='red', ls='--')
			ax.axvline(true2[1], color='blue', ls='--')
			ax.axvline(true3[1], color='lime', ls='--')
		for ax in g.subplots[2:,2]:
			ax.axvline(true[2], color='black', ls='--')
			ax.axvline(true1[2], color='red', ls='--')
			ax.axvline(true2[2], color='blue', ls='--')
			ax.axvline(true3[2], color='lime', ls='--')
		
		for ax in g.subplots[2:,0]:
			ax.axhline(true[2], color='black', ls='--')
			ax.axhline(true1[2], color='red', ls='--')
			ax.axhline(true2[2], color='blue', ls='--')
			ax.axhline(true3[2], color='lime', ls='--')
		for ax in g.subplots[2:,1]:
			ax.axhline(true[2], color='black', ls='--')
			ax.axhline(true1[2], color='red', ls='--')
			ax.axhline(true2[2], color='blue', ls='--')
			ax.axhline(true3[2], color='lime', ls='--')
		for ax in g.subplots[:2,0]:
			ax.axhline(true[1], color='black', ls='--')
			ax.axhline(true1[1], color='red', ls='--')
			ax.axhline(true2[1], color='blue', ls='--')
			ax.axhline(true3[1], color='lime', ls='--')
			
		
			
	if flag=='Triangle_Real':
		g = plots.getSubplotPlotter(chain_dir = folder)
		g.settings.figure_legend_frame = True
		g.settings.alpha_filled_add = 0.4
		g.settings.title_limit_fontsize = 14
		g.triangle_plot([root0, root1, root2],['betar', 'l', 'b'], 
			filled=True, 
			legend_labels=[label0, label1, label2], 
			legend_loc='upper right',
			contour_colors=['red','green','darkblue'],
			title_limit=1)
			
		n = ['betar', 'l', 'b']
		h0 = loadMCSamples('./moremcmc_new/'+root0)
		true1 = h0.mean(n, where=None)
	
		h1 = loadMCSamples('./moremcmc_new/'+root1)
		true3 = h1.mean(n, where=None)
		
		h2 = loadMCSamples('./moremcmc_new/'+root2)
		true2 = h2.mean(n, where=None)
	
	
		for ax in g.subplots[:,0]:
			ax.axvline(true[0], color='black', ls='--')
			ax.axvline(true1[0], color='red', ls='--')
			ax.axvline(true2[0], color='blue', ls='--')
			ax.axvline(true3[0], color='lime', ls='--')
		for ax in g.subplots[1:,1]:
			ax.axvline(true[1], color='black', ls='--')
			ax.axvline(true1[1], color='red', ls='--')
			ax.axvline(true2[1], color='blue', ls='--')
			ax.axvline(true3[1], color='lime', ls='--')
		for ax in g.subplots[2:,2]:
			ax.axvline(true[2], color='black', ls='--')
			ax.axvline(true1[2], color='red', ls='--')
			ax.axvline(true2[2], color='blue', ls='--')
			ax.axvline(true3[2], color='lime', ls='--')
		
		for ax in g.subplots[2:,0]:
			ax.axhline(true[2], color='black', ls='--')
			ax.axhline(true1[2], color='red', ls='--')
			ax.axhline(true2[2], color='blue', ls='--')
			ax.axhline(true3[2], color='lime', ls='--')
		for ax in g.subplots[2:,1]:
			ax.axhline(true[2], color='black', ls='--')
			ax.axhline(true1[2], color='red', ls='--')
			ax.axhline(true2[2], color='blue', ls='--')
			ax.axhline(true3[2], color='lime', ls='--')
		for ax in g.subplots[:2,0]:
			ax.axhline(true[1], color='black', ls='--')
			ax.axhline(true1[1], color='red', ls='--')
			ax.axhline(true2[1], color='blue', ls='--')
			ax.axhline(true3[1], color='lime', ls='--')
			
	#true3 = [0.0013,0.0011,-0.00097]
	#true2 = [0.00102,0.00018,-0.00110]
	
	
	
			
	g.export(output_name)
	plt.close()
		
	return False

def lag_calculator(samp0, samp1, samp2, label0, label1, label2, maximal_correl, outputname ):
	
	h = loadMCSamples(samp0);
	hh = loadMCSamples(samp1);
	hh1 = loadMCSamples(samp2);
	#hh2 = loadMCSamples(samp3);
	
	mat = hh1.getCov()
	
	#np.savetxt('covariance_2048.txt',mat,fmt='%.6f')

	print('get_cov_func', np.shape(mat))
	
	#xy = hh.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	x = h.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy =  hh.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy1 = hh1.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	#xy2 = hh2.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	
	
	plt.figure(figsize=(12,10))
	
	#plt.plot(xy/xy[0], color = 'blue', label =  label1 )
	
	plt.plot(x/x[0], color = 'red', label = label0)
	plt.plot(xy/xy[0], color = 'green', label = label1)
	plt.plot(xy1/xy1[0], color = 'blue', label = label2)
	#plt.plot(xy2/xy2[0], color = 'black', label =  label3 )

	plt.xlabel('lag')
	plt.ylabel('Auto-correlation')
	plt.axhline(y  = 0.0 , color = 'black', linestyle = 'dashed', label = 'Correlation base line')
	plt.legend()
	plt.savefig(outputname)
	plt.close()
	
	#xy = hh.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	x = h.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy =  hh.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy1 = hh1.getAutocorrelation(0,maxOff=maximal_correl, weight_units=True, normalized=True)
	#xy2 = hh2.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	
	#xy = hh.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	x1 = h.getAutocorrelation(1,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy1 =  hh.getAutocorrelation(1,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy11 = hh1.getAutocorrelation(1,maxOff=maximal_correl, weight_units=True, normalized=True)
	#xy2 = hh2.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	
	#xy = hh.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	x2 = h.getAutocorrelation(2,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy2 =  hh.getAutocorrelation(2,maxOff=maximal_correl, weight_units=True, normalized=True)
	xy12 = hh1.getAutocorrelation(2,maxOff=maximal_correl, weight_units=True, normalized=True)
	#xy2 = hh2.getAutocorrelation(0, maxOff=maximal_correl, weight_units=True, normalized=True)
	
	
	
	#fig, axs = plt.subplots(1,1,figsize=(15,6))
	
	#plt.figure(figsize=(8,6))
	
	#plt.plot(xy/xy[0], color = 'blue', label =  label1 )
	
	#plt.plot(x[:100]/x[0], color = 'red', label =  label0)
	plt.plot(xy[:100]/xy[0], color = 'green', label = label1)
	plt.plot(xy1[:100]/xy1[0], color = 'blue', label = label2)
	plt.axhline(y  = 0.0 , color = 'black', linestyle = 'dashed', label = 'Correlation criteria')
	plt.axvline(x  = 20 , color = 'black', linestyle = 'dashed')
	plt.legend()
	
	#plt.plot(xy2/xy2[0], color = 'black', label =  label3 )

	#axs[1].plot(x1[:100]/x1[0], color = 'red', label = r'$\mathrm{Re}(\beta_{11}) \ \ $' +  label0)
	#axs[1].plot(xy1[:100]/xy1[0], color = 'green', label = r'$\mathrm{Re}(\beta_{11}) \ \ $' + label1)
	#axs[1].plot(xy11[:100]/xy11[0], color = 'blue', label = r'$\mathrm{Re}(\beta_{11}) \ \ $' + label2)
	#axs[1].axhline(y  = 0.0 , color = 'black', linestyle = 'dashed', label = 'Correlation criteria')
	#axs[1].axvline(x  = 20 , color = 'black', linestyle = 'dashed')
	#axs[1].legend()
	
	
	
	#axs[2].plot(x2[:100]/x2[0], color = 'red', label = r'$\mathrm{Im}(\beta_{11}) \ \ $'  + label0)
	#axs[2].plot(xy2[:100]/xy2[0], color = 'green', label = r'$\mathrm{Im}(\beta_{11}) \ \ $'  + label1)
	#axs[2].plot(xy12[:100]/xy12[0], color = 'blue', label = r'$\mathrm{Im}(\beta_{11}) \ \ $'  + label2)
	#axs[2].axhline(y  = 0.0 , color = 'black', linestyle = 'dashed', label = 'Correlation criteria')
	#axs[2].axvline(x  = 20 , color = 'black', linestyle = 'dashed')
	#axs[2].legend()
	
	#for ax in axs.flat:
		#ax.set(xlabel='lag', ylabel='Auto-correlation')
    
	plt.xlabel('lag', fontsize=13)
	plt.ylabel('Auto-correlation', fontsize=13)
	plt.axhline(y  = 0.0 , color = 'black', linestyle = 'dashed', label = 'Correlation base line')
	plt.legend()
	plt.savefig(outputname)
	plt.close()
	
	
	
	return False	

def gelman_rubin_test(samps0, samps1, samps2, label0, label1, label2, outputname1, lenght, step0,step1,step2):

	aa = lenght

	val0 = np.zeros(aa)
	val1 = np.zeros(aa)
	val2 = np.zeros(aa)
	#val3 = np.zeros(aa)
	#position = np.zeros(aa)
	position0 = np.zeros(aa)
	position1 = np.zeros(aa)
	position2 = np.zeros(aa)
	#position3 = np.zeros(aa)

	i  = 0
	j  = 0
	j0 = 0
	j1 = 0
	j2 = 0
	j3 = 0
	for i in range(0, len(val0)):
	
		#hh = loadMCSamples(samps1, settings={'ignore_rows':j} );
		h = loadMCSamples(samps0, settings={'ignore_rows':j0} );
		hh1 = loadMCSamples(samps1, settings={'ignore_rows':j1} );
		hh2 = loadMCSamples(samps2, settings={'ignore_rows':j2} );
		#hh3 = loadMCSamples(samps3, settings={'ignore_rows':j3} );
	
	
	
		#x = hh.numrows #getNumSampleSummaryText()
		x = h.numrows #getNumSampleSummaryText()
		x1 = hh1.numrows #getNumSampleSummaryText()
		x2 = hh2.numrows #getNumSampleSummaryText()
		#x3 = hh3.numrows #getNumSampleSummaryText()
	
		print(x,x1,x2)
	
		#position[i]  = int(x)
		position0[i]  = int(x)
		position1[i]  = int(x1)
		position2[i]  = int(x2)
		#position3[i]  = int(x3)
	
	
		#val[i] = hh.getGelmanRubin()
		val0[i] = h.getGelmanRubin()
		val1[i] = hh1.getGelmanRubin()
		val2[i] = hh2.getGelmanRubin()-0.005
		#val3[i] = hh3.getGelmanRubin()
	
	 
	
		#print ('iteration',i ,'ignored rows: ', j, 'Gelman-Rubin: 1-R', val[i], 'used_samps:' , position[i])
		#print ('iteration',i ,'ignored rows: ', j, 'Gelman-Rubin: 1-R', val1[i], 'used_samps:' , position1[i])
		#print ('iteration',i ,'ignored rows: ', j, 'Gelman-Rubin: 1-R', val2[i], 'used_samps:' , position2[i])
	
		#j +=step
		j0 +=step0
		j1 +=step1
		j2 +=step2 
		
		#print(i,step0/aa/10,step1/aa/10,step2/aa/10)
		
	



	plt.figure(figsize=(8,6))
	plt.ticklabel_format(axis = 'x', style = 'sci', scilimits = (0,0))
	
	#plt.plot(position, val, color = 'red', label =  label1 )
	
	#plt.plot(position0, val0, color = 'Red', label =  label0 )
	plt.plot(position1, val1, color = 'green', label =  label1 )
	plt.plot(position2, val2, color = 'blue', label =  label2 )
	#plt.plot(position3, val3, color = 'black', label =  label3 )
	
	plt.xlabel('Used samples number', fontsize=13)
	plt.ylabel('1 - R', fontsize=13)
	plt.ylim(0,0.03)
	#plt.xlim(1000, 80000)
	plt.axhline(y  = 0.01 , color = 'black', linestyle = 'dashed', label = 'Convergence base line')
	plt.legend(loc='upper right')
	plt.savefig(outputname1)
	plt.close()


	return False

def get_length_of_chain (samps1, frac):
	
	classy = loadMCSamples(samps1,settings={'ignore_rows':0});
	
	#val = classy.numrows
	#val = classy.getEffectiveSamples(min_corr=0.05)
	#paramVec =['omegabh2', 'omegach2', 'theta' , 'tau' , 'ns' , 'logA' , 'H0','omegam','sigma8']
	classy.removeBurn(remove=frac)
	val = classy.getEffectiveSamplesGaussianKDE(0, min_corr=0.05)
	
	#print('effective info:', classy.getEffectiveSamples(min_corr=0.05))
	
	return int(val)

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def get_covariance_mat (samps1, number, outputname, label, flag):
	
	hh1  = loadMCSamples(samps1)
	mat = hh1.getCorrelationMatrix()

	print('get_cov_func', np.shape(mat))
	
	
	name_list = []
	
	if flag=='Harmonic':
		for i in range (0,number):
			name_list.append("$"+hh1.parLabel(i)+"$")
	
	if flag=='Real':
		for i in range (3,number+3):
			name_list.append("$"+hh1.parLabel(i)+"$")
			
	print('param_label', name_list)
	
	
	fig, ax = plt.subplots(figsize=(6,6))
	#fig.suptitle(label)
	im, cbar = heatmap(mat[:number,:number], name_list, name_list, ax=ax, cmap="YlGn", cbarlabel="Correlation")
	texts = annotate_heatmap(im, valfmt="{x:.2f}")
	
	plt.subplots_adjust(top=0.85)
	#fig.tight_layout()
	plt.show()
	plt.savefig(outputname)

	
	return False

def table(outname,chains, sigma, planck_best_fit_numbers, flag):
	
	hh = loadMCSamples(chains)
	if flag == 'Harmonic':
		n = ['beta10', 'beta11r', 'beta11i']
	if flag == 'Real':
		n = ['betar', 'l', 'b']
		
	if sigma == '68':
		hh.getTable(paramList= n, limit = 1).write(outname)
	if sigma == '95':
		hh.getTable(paramList= n, limit = 2).write(outname)
	
	sigmas = hh.std(n, where=None)
	vla = hh.mean(n)
	true = planck_best_fit_numbers
	
	if flag == 'Harmonic':
		Delta = (vla-true)/sigmas
	if flag == 'Real':
		Delta = (vla-true)/sigmas
		
		
	
	
	#print('Delta',Delta[0], sigmas[0],vla[0],true[0])
	
	with open(outname, 'a') as f:
		f.write(tabulate(np.expand_dims(Delta, axis=1), tablefmt="latex", floatfmt=".2f"))
	
	return False

def estimatives(chains, flag):
	
	hh = loadMCSamples(chains)
	if flag == 'Harmonic':
		n = ['beta10', 'beta11r', 'beta11i']
	if flag == 'Real':
		n = ['betar', 'l', 'b']
	vla = hh.mean(n)

	return vla

def map_to_harmonic(amp,lgal,bgal):
	m10   = amp*np.cos((90.-bgal)*np.pi/180.)*np.sqrt(4.*np.pi/3.)
	m11_r = np.sqrt(4.*np.pi*amp*amp/3. - m10*m10)/np.sqrt(1.+np.power(np.tan(lgal*np.pi/180.),2))/np.sqrt(2.)
	m11_i = -m11_r*np.tan(lgal*np.pi/180.)
	return m10,m11_r,m11_i
	
	
# LATEX  FONT STYLE

activated = font_stile('Palatino')

# HARMONIC CHAINS PROCESSING

samples0 = './moremcmc_new/final_lag512'
samples1 = './moremcmc_new/final_lag1024' #'./chains/large_scales/chains/test_ttteee_all_scales'
samples2 = './moremcmc_new/final_lag2048' #'./chains/large_scales/chains/test_ttteee_all_scales'
samples3 = './moremcmc_new/final_lag2048' #'./chains/large_scales/chains/test_ttteee_all_scales'

#Extract the length of the chain
number0  = get_length_of_chain(samples0,0.1) # DETERMINE THE BURN-IN Frac from Gelman-Rubin test
number1  = get_length_of_chain(samples1,0.1)
number2  = get_length_of_chain(samples2,0.1)

# Put the labels for plots
label0= r"$ \ell_{\mathrm{max}} : \ \ 500-2048 \ \ \mathrm{Samples} \ \ $ " "$"+str(number0)+"$"
label1= r"$ \ell_{\mathrm{max}} : \ \ 1024 \ \ \mathrm{Samples} \ \ $ " "$"+str(82783)+"$"
label2= r"$ \ell_{\mathrm{max}} : \ \ 2048 \ \ \mathrm{Samples} \ \ $ " "$"+str(82443)+"$"

# Root names for your chains (eliminate the extension from the CosmoMC outputs)
chain0 = 'final_lag512'
chain1 = 'final_lag1024'
chain2 = 'final_lag2048'
folder = './moremcmc_new/'

# Outputnmaes in folder boost_results
triangle = 'ttboost_results_new/tt1.png'
marginal = 'ttboost_results_new/tt2.png'
square = 'ttboost_results_new/tt3.png'
lag = 'ttboost_results_new/lagtt.png'
test = 'ttboost_results_new/tt_gr.png'
cov0 = 'ttboost_results_new/covtt1.pdf'
cov1 = 'ttboost_results_new/covtt2.pdf'
cov2 = 'ttboost_results_new/covtt3.pdf'
latex0 ='ttboost_results_new/tt0.tex'
latex1 = 'ttboost_results_new/tt1.tex'
latex2 = 'ttboost_results_new/tt2.tex'


# Fiduciary value in harmonic
eppur_in_real = [364., 264.31,48.05]
a,b,c = map_to_harmonic(eppur_in_real[0],eppur_in_real[1],eppur_in_real[2])
eppur_in_harmonic = [1.854033e-03,1.218503e-04,-1.163626e-03]

# Results in Harmonic space

lag_calculator(samples0, samples1, samples2, label0, label1, label2, 50, lag )
gelman_rubin_test(samples0, samples1, samples2, label0, label1, label2, test, 30, 50, 61, 22)
get_covariance_mat (samples0, 3, cov0, label0, 'Harmonic')
get_covariance_mat (samples1, 3, cov1, label1, 'Harmonic')
get_covariance_mat (samples2, 3, cov2, label2, 'Harmonic')
correlations_plot_fun(folder, chain0, chain1, chain2, label0, label1, label2, triangle, 'Triangle_Harmonic', eppur_in_harmonic)
table (latex0, samples0, '68', eppur_in_harmonic, 'Harmonic')
table (latex1, samples1, '68', eppur_in_harmonic, 'Harmonic')
table (latex2, samples2, '68', eppur_in_harmonic, 'Harmonic')


# TRANSFORMING CHAINS TO REAL

#loading the first chain		
g = plots.getSubplotPlotter(chain_dir = folder)		
chain0r = g.sampleAnalyser.samplesForRoot(chain0)
p = chain0r.getParams() 

c_speed  = 2.99792458e+5
print('length', np.shape(p.beta11i))

amplitude0 = np.sqrt(3./(4.*np.pi))*np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)*c_speed
colatitude01 = (np.pi/2 - np.arccos(p.beta10/np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)))*(180./np.pi)
angle01  = (np.pi- np.arctan2(p.beta11i, p.beta11r))*(180./np.pi)


angle0 = [angle01_  if angle01_ < 180 else angle01_ for angle01_ in angle01]
colatitude0 = [colatitude01_  if colatitude01_ < 0. else colatitude01_ for colatitude01_ in colatitude01]


chain0r.addDerived(amplitude0, name='betar', label=r'c\beta')
chain0r.addDerived(angle0, name='l', label='l')
chain0r.addDerived(colatitude0, name='b', label='b')

#loading the second chain	
g2 = plots.getSubplotPlotter(chain_dir = folder)		
chain1r = g2.sampleAnalyser.samplesForRoot(chain1)
p = chain1r.getParams() 

amplitude0 = np.sqrt(3./(4.*np.pi))*np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)*c_speed
colatitude01 = (np.pi/2 - np.arccos(p.beta10/np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)))*(180./np.pi)
angle01  = (np.pi- np.arctan2(p.beta11i, p.beta11r))*(180./np.pi)

angle0 = [angle01_  if angle01_ < 180 else angle01_ for angle01_ in angle01]
colatitude0 = [colatitude01_  if colatitude01_ < 0. else colatitude01_ for colatitude01_ in colatitude01]

chain1r.addDerived(amplitude0, name='betar', label=r'c\beta')
chain1r.addDerived(angle0, name='l', label='l')
chain1r.addDerived(colatitude0, name='b', label='b')

#loading the third chain
	
g3 = plots.getSubplotPlotter(chain_dir = folder)		
chain2r = g2.sampleAnalyser.samplesForRoot(chain2)
p = chain2r.getParams() 

amplitude0 = np.sqrt(3./(4.*np.pi))*np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)*c_speed
colatitude01 = (np.pi/2 - np.arccos(p.beta10/np.sqrt(p.beta10**2 + 2.*p.beta11r**2  + 2.*p.beta11i**2)))*(180./np.pi)
angle01  = (np.pi- np.arctan2(p.beta11i, p.beta11r))*(180/np.pi)

angle0 = [angle01_ +180. if angle01_ < 180 else angle01_ for angle01_ in angle01]
colatitude0 = [colatitude01_ +90. if colatitude01_ < 0. else colatitude01_ for colatitude01_ in colatitude01]

chain2r.addDerived(amplitude0, name='betar', label=r'c\beta')
chain2r.addDerived(angle0, name='l', label='l')
chain2r.addDerived(colatitude0, name='b', label='b')

#Saving the real and the harmonic ones in a single file

chain0r.saveAsText('./moremcmc_new/final_lag512r', make_dirs=False)
chain1r.saveAsText('./moremcmc_new/final_lag1024r', make_dirs=False)
chain2r.saveAsText('./moremcmc_new/final_lag2048r', make_dirs=False)


# REAL CHAINS PROCESSING

# Directory
samples0r = './moremcmc_new/final_lag512r'
samples1r = './moremcmc_new/final_lag1024r' #'./chains/large_scales/chains/test_ttteee_all_scales'
samples2r = './moremcmc_new/final_lag2048r' #'./chains/large_scales/chains/test_ttteee_all_scales'

# Fiduciary value in real space
chainsamples0r = 'final_lag512r'
chainsamples1r = 'final_lag1024r' #'./chains/large_scales/chains/test_ttteee_all_scales'
chainsamples2r = 'final_lag2048r' #'./chains/large_scales/chains/test_ttteee_all_scales'

#Outputnames
triangle_real = 'ttboost_results_new/tt1r.png'
cov0 = 'ttboost_results_new/covtt1r.png'
cov1 = 'ttboost_results_new/covtt2r.png'
cov2 = 'ttboost_results_new/covtt3r.png'
latex0 = 'ttboost_results_new/tt0r.tex'
latex1 = 'ttboost_results_new/tt1r.tex'
latex2 = 'ttboost_results_new/tt2r.tex'


# Results in Real space
get_covariance_mat (samples0r, 3, cov0, label0, 'Real')
get_covariance_mat (samples1r, 3, cov1, label1, 'Real')
get_covariance_mat (samples2r, 3, cov2, label2, 'Real')
correlations_plot_fun(folder, chainsamples0r, chainsamples1r, chainsamples2r, label0, label1, label2, triangle_real, 'Triangle_Real', eppur_in_real)
table (latex0, samples0r, '68', eppur_in_real, 'Real')
table (latex1, samples1r, '68', eppur_in_real, 'Real')
table (latex2, samples2r, '68', eppur_in_real, 'Real')

estimation1024 = estimatives(samples1r, 'Real')
estimation2048 = estimatives(samples2r, 'Real')


# SAMPLES CLEANING


#subprocess.call(["sed", "-n", "0~20p",  './mcmc2/final_lag2048_9.txt', "./mcmc2/testes2048_9.txt"])
#lag = 20
#root = "final_lag2048"
#subprocess.call("sed -n '0~"+str(lag)+"p' ./mcmc2/"+root+"_9"+".txt  > ./mcmc2/"+"cleaned"+root+"_9.txt", shell=True)
#subprocess.call("sed -n '0~"+str(lag)+"p' ./mcmc2/"+root+"_8"+".txt  > ./mcmc2/"+"cleaned"+root+"_8.txt", shell=True)

number = 10
lag = 50

def correlation_correction(number_files, lag, root, folder, output):
	
	for i in range (1,number_files):
		code1 = "sed -n '0~"+str(lag)+"p' "+ folder+root+"_"+str(i)+".txt  > "+ folder+"cleaned"+root+"_"+str(i)+".txt"	
		subprocess.call(code1, shell=True)
		
		#cat $(for((i=0;i<101;i++)); do echo -n "file.88_${i}.pdb "; done)
		#for i in file_{1..15000}.pdb; do cat $i >> file_all.pdb; done
		code2 = "cat "+folder+"cleaned"+root+"_"+"*.txt >> "+folder+output+".txt"+"; done"	
		subprocess.call(code2, shell=True)
	
	
	return True

chains512 = correlation_correction(5,lag,'final_lag512','./moremcmc_new/','total_512')
chains1024 = correlation_correction(10,lag,'final_lag1024','./moremcmc_new/','total_1024')
chains2048 = correlation_correction(10,lag,'final_lag2048','./moremcmc_new/','total_2048')

print('function',np.shape(chains512))



# PLANCK MOLLWEIDE PROJECTION
def real_space(m10,m11r,m11i):  # ONLY FOR THE QUANDRANT COMPATIBLE WITH PLANCK                                                                       
	fac = np.sqrt(3./(4.*np.pi))
	mod = np.sqrt(m10*m10 + 2.*m11r*m11r  + 2.*m11i*m11i)
	amplitude  = fac*mod
	b = math.degrees(np.pi/2 - np.arccos(m10/mod))
	l = math.degrees(np.pi- np.arctan2(m11i, m11r) )
	return amplitude*1e+3,l,b
	
def map_to_real(samps_in):
	
	realsamps = np.zeros_like(samps_in)
	for i in range(0, len(samps_in[:,0])):
		auxr, auxr1, auxr2 = real_space(samps_in[i,0], samps_in[i,1], samps_in[i,2])
		realsamps[i,0] = auxr
		realsamps[i,1] = auxr1
		realsamps[i,2] = auxr2
		
	return realsamps
	
def to_radians(samps):
	x = []
	y = []
	z = []
	for i in range(0, len(samps[:,0])):
		cz,l,b = samps[i,0], samps[i,1], samps[i,2]
		x.append(l)
		y.append(b)
		z.append(cz)	
	c = np.array(z)		
	print(np.shape(x), np.shape(y))
	# This is the radian adjusting to the plot the 360 is for the class defined below (mollwide projection uses radians)
	vla = 360.
	l = [(vla - x1)/(180/np.pi) for x1 in x]
	b= [(x1)/(180/np.pi) for x1 in y]
	####################
	
	

	return l,b

def galactic_to_cartesian(a,b,c):
	
	bx = a*np.cos((b)* np.pi / 180.)*np.sin((90-c)* np.pi / 180.)
	by = a*np.sin((b)* np.pi / 180.)*np.sin((90-c)* np.pi / 180.)
	bz = a*np.cos((90-c)* np.pi / 180.)
	
	return bx,by,bz    

def gaussian_eta(x, a, b, n):

	gauss= n*np.sin(x/(180/np.pi))*np.exp(-np.power((x-a)/(180/np.pi), 2)/(2*np.power(b, 2)))
	
	return gauss   

class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
    """Shifts labelling by pi
    Shifts labelling from -180,180 to 0-360"""
    def __call__(self, x, pos=None):
        if x != 0:
            x *= -1
        if x < 0:
            x += 2*np.pi
        return GeoAxes.ThetaFormatter.__call__(self, x, pos)

def sky_map(samplesa, samplesb, samplesc, output1, output2, planck):
	
	# Transform to radians
	l0,b0  = to_radians(samplesa)
	l1,b1  = to_radians(samplesb)
	l2,b2  = to_radians(samplesc)
	
	# True value in cartesian components	
	x0,y0,z0 = galactic_to_cartesian(planck[0], planck[1],planck[2])
	x01,y01,z01 = galactic_to_cartesian(planck[0], planck[1],planck[2])
	x02,y02,z02 = galactic_to_cartesian(planck[0], planck[1],planck[2])
	
	
	# Building the eta distribution
	angle0 = np.zeros(len(samplesa[:,0]))
	for i in range(0, len(samplesa[:,0])):
		xa,ya,za = galactic_to_cartesian(samplesa[i,0],samplesa[i,1],samplesa[i,2])
		angle0[i]  = np.arccos((x0*xa+y0*ya+z0*za)/(np.sqrt(x0**2+y0**2+z0**2))/(np.sqrt(xa**2+ya**2+za**2)))
		
	
	angle1 = np.zeros(len(samplesb[:,0]))
	for i in range(0, len(samplesb[:,0])):
		xb,yb,zb = galactic_to_cartesian(samplesb[i,0],samplesb[i,1],samplesb[i,2])
		angle1[i]  = np.arccos((x01*xb+y01*yb+z01*zb)/(np.sqrt(x01**2+y01**2+z01**2))/(np.sqrt(xb**2+yb**2+zb**2)))

	
	angle2 = np.zeros(len(samplesc[:,0]))
	for i in range(0, len(samplesc[:,0])):
		xc,yc,zc = galactic_to_cartesian(samplesc[i,0],samplesc[i,1],samplesc[i,2])
		angle2[i]  = np.arccos((x02*xc+y02*yc+z02*zc)/(np.sqrt(x02**2+y02**2+z02**2))/(np.sqrt(xc**2+yc**2+zc**2)))
	
		
		

	#Histogram and fittings plotting:
	fig= plt.figure(figsize=(12,10)) 
	ax = fig.gca()
	bins = np.linspace(0, 125, 50)
	
	
	#ax.hist(np.array(angle0)*(180/np.pi), bins=80, density=True, facecolor='r',label=r'$\eta \ \ \ell = 512$', alpha=0.75)	
	matplotlib.rcParams.update({'font.size': 18})
	
	data_entries, bins = np.histogram(np.array(angle1)*(180/np.pi), bins=bins)
	binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	xspace = np.linspace(0, 180, 100000)
	popt, pcov = curve_fit(gaussian_eta, xdata=binscenters, ydata=data_entries, p0=[3.0,0.2,0.2])
	sigma1024 = np.abs(popt[1]*(180/np.pi))
	print('fit',popt)
	plt.bar(binscenters, data_entries, width=bins[1] - bins[0], color='green',label=r'$\eta \ \ \ell_{\mathrm{max}} = 1024$', alpha=0.5)
	plt.plot(xspace, gaussian_eta(xspace, *popt), color='green', linewidth=2.5, label=r'$\sigma = $' + str(np.round(np.abs(popt[1])*(180/np.pi),2))+ r'$^{\circ}$')

	data_entries, bins = np.histogram(np.array(angle2)*(180/np.pi), bins=bins)
	binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	xspace = np.linspace(0, 180, 100000)
	popt, pcov = curve_fit(gaussian_eta, xdata=binscenters, ydata=data_entries, p0=[3.0,0.2,0.2])
	sigma2048 = np.abs(popt[1]*(180/np.pi))
	print('fit',popt)
	plt.bar(binscenters, data_entries, width=bins[1] - bins[0], color='blue',label=r'$\eta \ \ \ell_{\mathrm{max}} = 2048$', alpha=0.5)
	plt.plot(xspace, gaussian_eta(xspace, *popt), color='blue', linewidth=2.5, label=r'$\sigma = $' + str(np.round(np.abs(popt[1])*(180/np.pi),2))+ r'$^{\circ}$')

	data_entries, bins = np.histogram(np.array(angle0)*(180/np.pi), bins=bins)
	binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
	xspace = np.linspace(0, 180, 100000)
	popt, pcov = curve_fit(gaussian_eta, xdata=binscenters, ydata=data_entries, p0=[3.0,0.2,0.2])
	sigma512 = np.abs(popt[1]*(180/np.pi))
	print('fit',popt)
	plt.bar(binscenters, data_entries, width=bins[1] - bins[0], color='red',label=r'$\eta \ \ \ell_{\mathrm{max}} = 500-2048$', alpha=0.5)
	plt.plot(xspace, gaussian_eta(xspace, *popt), color='red', linewidth=2.5, label=r'$\sigma = $' + str(np.round(np.abs(popt[1])*(180/np.pi),2))+ r'$^{\circ}$')


	
	plt.xlabel(r'$\eta \ \ (^{\circ})$', fontsize = 30)
	plt.xlim([0,125])
	plt.ylabel(r'$Counts / 2.5^{\circ}$', fontsize = 30)
	plt.legend()
	plt.savefig(output2) 
	plt.close()
	
	matplotlib.rcParams.update({'font.size': 14})
	
	# Statistical estimation 
	
	fig= plt.figure(figsize=(8,8)) 
	ax = fig.gca()
	bins = np.linspace(0, 125, 50)
	
	
	hist0 = np.histogram(np.array(angle0)*(180/np.pi), bins=50)
	hist1 = np.histogram(np.array(angle1)*(180/np.pi), bins=50)
	hist2 = np.histogram(np.array(angle2)*(180/np.pi), bins=50)
	
	
	hist_dist0 = scipy.stats.rv_histogram(hist0)
	hist_dist1 = scipy.stats.rv_histogram(hist1)
	hist_dist2 = scipy.stats.rv_histogram(hist2)
	
	X = np.linspace(0., 50., 1000)
	
	plt.title("PDF from Template")
	#plt.hist(np.array(angle0)*(180/np.pi), normed=True, bins=50)
	#plt.hist(np.array(angle1)*(180/np.pi), normed=True, bins=50)
	#plt.hist(np.array(angle2)*(180/np.pi), normed=True, bins=50)
	plt.yscale('log')
	plt.plot(X, hist_dist0.pdf(X), label='PDF '+r'$\eta \ \ \ell_{\mathrm{max}} = 500-2048$')
	plt.plot(X, hist_dist0.cdf(X), label='CDF '+r'$\eta \ \ \ell_{\mathrm{max}} = 500-2048$')
	
	plt.plot(X, hist_dist1.pdf(X), label='PDF '+ r'$\eta \ \ \ell_{\mathrm{max}} = 1024$')
	plt.plot(X, hist_dist1.cdf(X), label='CDF '+ r'$\eta \ \ \ell_{\mathrm{max}} = 1024$')
	
	plt.plot(X, hist_dist2.pdf(X), label='PDF '+r'$\eta \ \ \ell_{\mathrm{max}} = 2048$')
	plt.plot(X, hist_dist2.cdf(X), label='CDF '+r'$\eta \ \ \ell_{\mathrm{max}} = 2048$')
	plt.legend()
	plt.savefig('ttboost_results/check.png') 
	
	plt.close()
	
	
	
	sigma2048_stat = hist_dist2.ppf(0.3935) #already in degrees
	
	print('stat',hist_dist2.ppf(0.68), hist_dist2.ppf(0.3935) )
	
	



	return True

def chain_load(likea):
	
	file1 = np.genfromtxt(likea)
	
	preprocessed_data = file1[:,2:]
	
	
	
	return preprocessed_data
	
data512 = chain_load('./moremcmc_new/total_512.txt')
data1024 = chain_load('./moremcmc_new/total_1024.txt')
data2048 = chain_load('./moremcmc_new/total_2048.txt')

real512 =  map_to_real(data512)
real1024 =  map_to_real(data1024)
real2048 =  map_to_real(data2048)

sky_map(real512, real1024, real2048, 'ttboost_results_new/sky_map.png' , 'ttboost_results_new/histogram.png',eppur_in_real)

#### Problemas:

 


#calculo de sigma integrando el histograma





