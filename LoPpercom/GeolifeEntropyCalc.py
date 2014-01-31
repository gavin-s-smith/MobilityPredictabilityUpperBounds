'''
Created on 10 Jul 2013

@author: Romain Wieser
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




'''

from __future__ import division
from datetime import timedelta
from GeolifeSymbolisation import get_geolife_data
from mlabwrap import mlab # @UnresolvedImport This is the import for mlabwrap
import numpy as np
import pylab as pl
import time
from Utils import ensure_dir
from GenericLoP import empiricalEntropyRate
import GeolifeSymbolisation

def parse_timedelta(time_str):
    """
    Helper function to convert human readable time offsets to Python timedelta objects.
    This is simply to enable readability in the parameter specification.
    
    :param time_str: A time offset string in the form 'h:mm:ss'
    :type time_str: str
    """
    split = time_str.split(':')
    return timedelta(hours = int(split[0]), minutes = int(split[1]), seconds = int(split[2]))

####################################################################################################
######################              Some hard-coded parameters              ########################
####################################################################################################

listSpatialRes = [1000,5000,10000,50000,100000,1000000,3000000,10000000,50000000]
listRealSpatialRes = [618,2474,9896,39586,158347,633388,2533555,10134220,40536880]
#listSpatialRes = [5000]
#listRealSpatialRes = [2474]
dict_SR = dict(zip(listSpatialRes,listRealSpatialRes))
listTemporalRes = ['0:05:00','0:10:00','0:15:00','0:30:00','0:45:00','1:00:00']
#listTemporalRes = ['0:05:00']
listTemporalRes = [parse_timedelta(temporalRes) for temporalRes in listTemporalRes]
listTemporalResSecond = [temporalRes.total_seconds() for temporalRes in listTemporalRes]

#####################################################################################################

def save_results( file_name, LoP, DL_RL):
    """
    Writes the upper bounds computed for the upper limit of predictability to a CSV file.
    
    :param file_name: Base filename to write (minus DL/RL type and .csv)
    :type file_name: str
    :param LoP: 2D array of upper bound values
    :type LoP: 2D array of reals
    :param DL_RL: String denoting original method (DL) or refined method (RL)
    :type DL_RL: str
    """
    #Save as a csv file :
    f = file(file_name +DL_RL+ ".csv", 'w')
    np.savetxt(f, LoP,fmt ="%.5f")
    f.close()



def run( group = "All",scale = None, output_dir = './ResultsLoP_replication/final_graphs', bulk_build_preprocessing = False):
    """
    Generates a single heatmap for a given list of Geolife ids, for a given method of computing the upper bound on
    the upper limit of predictability.
    
    :param group: ["id_str",[list of ids in the geolife dataset]]
    :type group: Nested list
    :param scale: [min_z, max_z, step]  Set the scale of the heatmap z
    :type scale: Float array
    """
    t = time.time()
    
    #Group setting
    if(group == "All"):
        suffix = "All"
        persons = "All"
    else:
        suffix = "Grp{}".format(group[0])
        persons = group[1]
    
    
    if not output_dir[-1] == '/':
        output_dir = output_dir + '/'
        
    file_name = "{}Heatmap_{}".format(output_dir,suffix)
    
    ensure_dir(file_name)
    
    print "Calculing the LoP for {}".format(suffix)
    
    if bulk_build_preprocessing:
        # will attempt to bulk build the cache using multiple CPU cores
        # will skip caches if already built.
        # if this option is not specified and a cache does not exist
        # it will be built when required, using a single CPU core.
        GeolifeSymbolisation.bulk_build_resolution_cache(listSpatialRes, listTemporalRes)
    
    mlab.openPool()
    failed_ids = set()
    LoP_RL = []
    LoP_DL = []
    LoP_failed_ct = []
    passed_norm_test = []
    for spatialRes in listSpatialRes:
        LoP_RL.append([])
        LoP_DL.append([])
        LoP_failed_ct.append([])
        passed_norm_test.append([])
        for temporalRes in listTemporalRes:
            
            #Compute data

            #---------------------------------------------
            #Load data from an existing preproc database, this will have been created
            # earlier if it did not exist.    
            data, person_ids = get_geolife_data(spatialRes, temporalRes,persons)
            #---------------------------------------------
            
            # Sanity check on loading
            for person in data:
                if len(person) == 0:
                    raise Exception("One or more person's trajectory was not loaded/created correctly.")
            # End sanity check
            
            S_RL, N_RL = empiricalEntropyRate(data,'RL')
            S_DL, N_DL = empiricalEntropyRate(data,'DL')
                    
            #Save the average:

            tmpG_RL = list(mlab.ParLoP(S_RL, N_RL)[0])
            tmpG_DL = list(mlab.ParLoP(S_DL, N_DL)[0])
            
            #-88 real fail in solve
            #-99 known fail in solve when S > log2(N)
            # See the Matlab script (ParLoP.m) for more details
            
            if (np.asarray(tmpG_RL)==-88).any():
                raise Exception("ERROR: (RL) Matlab failed the solve, but the entropy was in the correct range. Therefore an unknown error has occured.")
            
            if (np.asarray(tmpG_DL)==-88).any():
                raise Exception("ERROR: (DL) Matlab failed the solve, but the entropy was in the correct range. Therefore an unknown error has occured.")
            
            
            # Replace known solve fails. These are the cases when an entropy is found that is to high. 
            # This means the LZ entropy rate estimate is wrong (the estimator has failed to converge)
            # There is no way to correct this, without collecting more data from the individual.
            # While excluding the individual is not ideal it is better than including a value that is 
            # *known* to be erroneous. Therefore we discard the individual. 
            tmpG_RL = np.asarray(tmpG_RL)
            tmpG_DL = np.asarray(tmpG_DL)
            
            tmpG_RL_known_fail_mask = tmpG_RL < -1
            tmpG_DL_known_fail_mask = tmpG_DL < -1
            

            # To be comparable we must arrive at a consistent set of individuals from which to compare both 
            # methods. 
            tmpG_known_fail_mask = np.asarray(tmpG_RL_known_fail_mask) | np.asarray(tmpG_DL_known_fail_mask)
            
            #print tmpG_known_fail_mask
            
            failed_ct = len(tmpG_RL[tmpG_known_fail_mask])
            
            for p in np.asarray(person_ids)[tmpG_known_fail_mask]:
                failed_ids.add(p)
            

            # Filter out known solve fails.
            tmpG_RL = list(np.asarray(tmpG_RL)[~tmpG_known_fail_mask])
            tmpG_DL = list(np.asarray(tmpG_DL)[~tmpG_known_fail_mask])
            
            if not len(tmpG_RL) == len(tmpG_DL):
                raise Exception("SHOULD NOT OCCUR 5g4dfg65")
            
            if (np.asarray(tmpG_RL) < 0).any():
                raise Exception("ERROR. lsdkfal")
            
                
            LoP_RL[-1].append(np.average(tmpG_RL))
            LoP_DL[-1].append(np.average(tmpG_DL))
            LoP_failed_ct[-1].append( failed_ct )
                

            
    mlab.closePool()
    

    
    save_results( file_name, LoP_RL, 'RL')
    save_results( file_name, LoP_DL, 'DL')
    
    f2 = file(file_name + "_failed_ct.csv", 'w')
    print 'failed_ids = {}.'.format( failed_ids )
    
    np.savetxt(f2, LoP_failed_ct,fmt ="%.5f")
    f2.close()
    
    print "Done in {} seconds".format(time.time() - t)
    
    

if __name__ == '__main__':
    """
    Replicates the results from the PERCOM 2014 paper:
    A Refined Limit on the Predictability of Human Mobility
    by Gavin Smith, Romain Wieser, James Goulding and Duncan Barrack
    
    """
    # Parameters
    output_dir = './ResultsLoP_replication/final_graphs'
    
    useable_ids =  ["PERCOM",[0, 1, 2, 3, 4, 5, 7, 9, 12, 13, 14, 15, 16, 17, 22, 24, 153, 28, 30, 35, 36, 38, 39, 40, 43, 44, 50, 179, 52, 55, 68, 71, 82, 84, 85, 92, 96, 101, 104, 167, 119, 126]]
    
    # Warning: If bulk_build_preprocessing is True, then the ONLY option for useable_ids is the PERCOM paper set.
    # This is because the bulk build code only caches the PERCOM paper IDs (this is hardcoded into the method).
    # Using a different set of useable_ids along with bulk_build_preprocessing=True WILL result in an exception being thrown
    # due to the IDs subsequently not being in the cache. If you would like to bulk build with a different ID set, change the
    # personIds list hardcoded in GeolifeSymbolisation.py::buildPreprocessingTable
    run(useable_ids, scale = [0.05,1,0.05], output_dir = output_dir, bulk_build_preprocessing=True) # scale controls the output headmap colour axis
    
    