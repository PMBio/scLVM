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

#plot formatting scripts and tools
#plot settings
import pylab as PL

plparams = {'backend': 'pdf',
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 13,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': False}
PL.rcParams.update(plparams)


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
