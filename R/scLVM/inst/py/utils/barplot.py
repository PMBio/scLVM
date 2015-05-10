# Copyright(c) 2014, The scLVM developers (Forian Buettner, Paolo Francesco Casale, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import scipy as SP
import sys
import os
sys.path.append('../')
import pdb
import glob
import cPickle
import h5py
import matplotlib as mpl
import pylab as PL
import scipy.stats as ST
import matplotlib.gridspec as gridspec
from plot_format import *

def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

def var_plot(var,H2,fields,normalize=False,plot_density=True,plot_element_count=False,V_range=SP.linspace(0,1,11),filename=None,figsize=[6,5]):

    """
    var:    GxC matrix, each row contains the variance components of gene g
    H2:     G vector used to bin the genes (in this example it is the total heritability)
    fields: array with field descriptors[index,legend,color]; see default.py for an illustration
    normalize: normalize variance in each bin to 1 (False)
    plot_element_count: plot number of elements in bin (True)
    plot_density: plot the gene density
    V_range: range of variance bins
    filename: figure.pdf for export (None)
    """

    #open figure, etc.
    PL.figure(figsize=figsize)
    if plot_density: 
        PL.axes([0.1,0.16,0.8,0.72])

    # plot thin horizonatal lines, every second bin of plot range
    PL.xlim(0,1)
    for y in V_range[::2]:
        PL.plot([0,1],[y,y],':',color='Gray',lw=0.1)
    
    _x = 0.01
    dx = 0.1
    wx = 0.08
    #plot headers for legend
    H = [None]*len(fields)
    #xticks / xticks_labels
    xticks = []
    xticks_labels = []
    for i in xrange(len(V_range)-1):
        #create label
        xticks.append(V_range[i:i+1].mean())
        xticks_labels.append('%d-%d%%' % tuple(V_range[i:i+2]*100))
        #variance bin
        v0 = V_range[i]
        v1 = V_range[i+1]
        #select genes matching
        Irel = (H2>=v0) & (H2<v1)
        N = Irel.sum()
        #average in bin 
        V = var[Irel,:].mean(0)
        if normalize:
            #fields that are plot?
            Iplot = SP.array([field[0] for field in fields])
            V[Iplot]/=V[Iplot].sum()

        #plot box plot for each plot descriptora
        _y = 0
        for i in xrange(len(fields)):
            field = fields[i]
            #field is [index,legend,color]
            #bar height
            bh   = V[field[0]]
            #plot
            H[i] =PL.bar(_x,bh,wx,bottom=_y,color=field[2])
            #increment
            _y += bh

        #number of genes in bin
        if plot_element_count:
            PL.text(_x+wx/2,_y+0.04,'(%d)' % (N),fontsize=8,horizontalalignment='center')
        #move to nex tbin
        _x +=dx

    #limit ylim to make sure we don't any issues..
    PL.ylim([0,1])
    # build the legend
    legends = [field[1] for field in fields]
    lh=PL.legend(H,legends,
              loc='upper center',
              bbox_to_anchor=(0.5, 1.21),
              ncol=5,
              prop={'size':10})
    lh.set_frame_on(False)

    # x axis tick labels
    plt = PL.gca()
    simpleaxis(plt)
    plt.set_xticks(xticks)
    plt.set_xticklabels(xticks_labels)
    PL.xticks(rotation=40)
    #fix yticks
    yticks = plt.get_yticks()
    yticks_labels = ['%d%%' % tick for tick in (yticks*100)]
    plt.set_yticklabels(yticks_labels)

    # remove top and right axis
    plt.spines["right"].set_visible(False)
    plt.spines["top"].set_visible(False)
    plt.xaxis.set_ticks_position('bottom')
    plt.yaxis.set_ticks_position('left')

    # axis labels
    PL.ylabel('Explained variance')
    PL.xlabel('Total non-technical variance')

    yticks = PL.yticks()[0][0:-1]
    PL.yticks(yticks)
    if plot_density:
        dens=ST.gaussian_kde(H2[~SP.isnan(H2)])
        x = SP.linspace(0,1,500)
        y = dens.evaluate(x)
        PL.axes([0.1,0.87,0.8,0.07])
        #PL.hist(h[~SP.isnan(h)])
        PL.plot(x,y,'k-')
        PL.xlim([0,1])
        PL.xticks([])
        ym = PL.ylim()[1]
        yticks = SP.array([0,0.5*ym,ym])
        PL.yticks(yticks, fontsize=8)
        PL.ylabel('density', fontsize=8)
        plt = PL.gca()
        simpleaxis(plt)

    # save figure
    if filename is not None:
        PL.savefig(filename, bbox_inches='tight')

        
