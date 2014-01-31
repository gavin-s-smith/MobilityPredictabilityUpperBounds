'''
Created on 22 Oct 2013

@author: Gavin Smith
@organization: Horizon Digital Economy Institute, The University of Nottingham.

@copyright: This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Dependencies:
* Matlab
* mlabwrap
* GPU library libEntropyCalc

'''

from __future__ import division
import numpy as np
from mlabwrap import mlab # @UnresolvedImport This is the import for mlabwrap
import libEntropyCalc as GPU_LZ_EC # @UnresolvedImport
import math

def get_N_DL(sym_list):
    """
    Computes the number of distinct locations in the trajectory
    :param sym_list: A list of location symbols
    :type sym_list: list
    """
    
    return len(np.unique(sym_list))

def get_N_RL(sym_list):
    """
    Compute a value denoting the maximum "number of reachable locations", ($N_{r}$), over all possible locations.
    
    From the paper:
    Formally $N_{r}$ is calculated from an empirical symbolic time series $\mathcal{T} = \{s_{1}, s_{2}, \ldots, s_{m}\}$, 
    with the set of all possible spatial locations being $\Omega$, as $N_{r} = \max_{x \in \Omega} | \{ s_{i+1} : s_i = x \} |$.
    
    :param sym_list: A list of location symbols
    :type sym_list: list
    """
    
    mapLocation = {}
    ct_point = 0
    for point in sym_list[:-1]:
        idNextPoint = sym_list[ct_point+1]
        try:
            mapLocation[point].add(idNextPoint)
        except KeyError:
            mapLocation[point] = set([idNextPoint])
        ct_point += 1
    N = 0
    for SetPoint in mapLocation.values():
        nb = len(SetPoint)
        if nb > N:
            N = nb
    return N
     
    
def empiricalEntropyRate(data,N_mode = "DL"):
    

    print "Computing empirical entropy rate..."
    
    empiricalEntropyRate = [] 
    N = []
    if(N_mode == "DL-RL"):
        N.append([])
        N.append([])
    
    for person in data :
        
        sym_list_orig = np.array(person)
        
        #Entropy resolution:
        #prepare data
        sym_list = sym_list_orig.astype(np.int64)
        output = np.zeros(len(sym_list_orig), dtype=np.int64)
        output = output.astype(np.int64)
        n = len(sym_list[sym_list >= 0])
        # Use EC lib
        GPU_LZ_EC.EC( sym_list, output  )
        #Calc the entropy :
        gpu_ent = math.pow(sum( [ output[i] / math.log(i+1,2) for i in range(1,n)] ) * (1.0/n),-1)
        #Append the entropy
        empiricalEntropyRate.append(gpu_ent)
        
        #N resolution:
        if(N_mode == "DL"):
        #------------------------=0 Distinct Location 0=-------------------------
            N.append(get_N_DL(sym_list))
        elif(N_mode == "RL"):
        #------------------------=0 Reachable Location 0=------------------------
            N.append(get_N_RL(sym_list))
        else:
            raise Exception( "Error: Unknown N_mode. Only DL or RL known, {} given.".format(N_mode) )

    print "S :", empiricalEntropyRate
    print "N : ", N
    
    return empiricalEntropyRate, N




def process_symbolic_data( data, standard_method = True, refined_method = False):
    """
    Given a list of trajectories (regularly sampled location integer symbols) returns
    the request upper bound(s) on the upper limit of predictability.
    
    Returns either a single value or two [standard_method, refined_method].
    
    :param data: List of trajectories
    :type data: List of List of int
    :param standard_method: True to calculate the upper bound on the upper limit of predictability using the original method by Song et. al. 
    :type standard_method: Boolean
    :param refined_method: True to calculate the upper bound on the upper limit of predictability using the refined method from our PERCOM paper.
    :type refined_method: Boolean
    """
    
    if refined_method:
        S_RL, N_RL = empiricalEntropyRate(data,'RL')
    
    if standard_method:
        S_DL, N_DL = empiricalEntropyRate(data,'DL')
                 
    mlab.openPool()
    
    if refined_method:
        tmpG_RL = list(mlab.ParLoP(S_RL, N_RL)[0])
    
    if standard_method:
        tmpG_DL = list(mlab.ParLoP(S_DL, N_DL)[0])
    
    mlab.closePool()
    
    if standard_method:
        print '\nStandard method: AVG: {} MIN: {} MAX: {}'.format(np.mean(np.asarray(tmpG_DL)),np.min(np.asarray(tmpG_DL)),np.max(np.asarray(tmpG_DL)))

    if refined_method:
        print '\nRefined method: AVG: {} MIN: {} MAX: {}'.format(np.mean(np.asarray(tmpG_RL)),np.min(np.asarray(tmpG_RL)),np.max(np.asarray(tmpG_RL)))

    if refined_method and standard_method:
        return tmpG_DL, tmpG_RL
    
    if refined_method:
        return tmpG_RL
    
    return tmpG_DL

