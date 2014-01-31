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
'''

from __future__ import division
import apsw
import numpy as np
import healpy as hp  # @UnresolvedImport
from datetime import datetime as dt
from Utils import ensure_dir
import os
from multiprocessing import Pool, cpu_count
from datetime import timedelta
fmt = '%Y-%m-%d %H:%M:%S'

#========
# PATHS
#========

main_geolifedb = '../DataGeolife/geolife.sqlite' # will be created if it doesn't exist
preprocessing_dir = '../DataGeolife/Preproc' # folder for semi-temporary databases for specific quantisations, will be created if they don't exist
geolife_zipfile_data_dir = '~/Downloads/Geolife Trajectories 1.3/Data' # folder for the original data. Used to build the above databases if they don't exist.



ensure_dir(main_geolifedb)
ensure_dir(preprocessing_dir)


if '~' in geolife_zipfile_data_dir:
    geolife_zipfile_data_dir = os.path.expanduser(geolife_zipfile_data_dir)

#========

def build_specific_cache( spatialRes_temporalRes_pair ):
    spatialRes = spatialRes_temporalRes_pair[0]
    temporalRes = spatialRes_temporalRes_pair[1]
    
    print 'Processing spatial res {}, temporal res {}'.format( spatialRes, temporalRes )
    
    if not os.path.exists( "{}/S{}T{}.sqlite".format(preprocessing_dir,spatialRes, temporalRes) ):
        buildPreprocessingTable(spatialRes,temporalRes,nest = True)
    

def bulk_build_resolution_cache(listSpatialRes, listTemporalRes, personsId = "All" ):
    """
    Bulk builds the spatial/temporal resolution cache files. This can be done on
    the fly but this can not be done in parallel. Since this operation takes a long
    time this method has been written to take advantage of multiple cores.
    
    :param listSpatialRes: A list of spatial resolutions required for the data.
    :type listSpatialRes: list of ints denoting meters
    :param listTemporalRes: A list of temporal resolutions required for the data.
    :type listTemporalRes: list of datetime.timedelta
    :param personsId: List of the person IDs for which the data should be fetched.
    :type personsId: List of ints
    """
    
    if not os.path.exists( main_geolifedb ):
        build_main_db()
    
    pairs = []
    for spatialRes in listSpatialRes:
        for temporalRes in listTemporalRes:
            pairs.append( (spatialRes, temporalRes) )
            
            #for debugging
            #build_specific_cache( pairs[-1] )
    #exit(-1)
    
    pool = Pool( processes = cpu_count() - 2 ) # leave some CPU for day to day tasks :-), 2 actual is one real CPU core on a Intel hyperthreaded system
    pool.map(build_specific_cache, pairs ) # the function "match( symbol_idx )" will now be called in parallel with the argument 0,1,... etc.
    pool.close()
    pool.join()
    



def get_geolife_data(spatialRes, temporalRes, personsId = "All" ):
    """
    Loads Geolife data for a given spatiotemporal resolution and a specific set of person IDs.
    Builds an SQLite table caching the quantisation from the original dataset if it does not exist.
    
    :param spatialRes: The spatial resolution required for the data.
    :type spatialRes: int denoting meters
    :param temporalRes: The temporal resolution required for the data.
    :type temporalRes: datetime.timedelta
    :param personsId: List of the person IDs for which the data should be fetched.
    :type personsId: List of ints
    """
    
    if not os.path.exists( "{}/S{}T{}.sqlite".format(preprocessing_dir,spatialRes, temporalRes) ):
        if not os.path.exists( main_geolifedb ):
            build_main_db()
        else:
            # ensure it has the table we need in it
            pass
            
        buildPreprocessingTable(spatialRes,temporalRes,nest = True, personsIds=personsId)
    
    return loadData(spatialRes, temporalRes, personsId )



#==============================================
#        Helper methods

#load preprocessing data  from the data base :
def loadData(spatialRes, temporalRes, personsId = "All", withDates = False):
    print "loading..."
    
    connection = apsw.Connection("{}/S{}T{}.sqlite".format(preprocessing_dir,spatialRes, temporalRes))
    
    curs1 = connection.cursor()
    curs2 = connection.cursor()
    curs3 = connection.cursor()
    
    data =[]
#    idPerson = []
    if(personsId == "All"):
        #All
        sql = "SELECT DISTINCT person FROM preproc GROUP BY person "
    #    #    WHERE person = 35
        personsId = map(lambda x: x[0],curs1.execute(sql) )
        
#    personsId = [49,137,143,149,151,177,106] #whith higher Q
#    personsId = [153,128,41,17,163,68,10] #with higher number of points
    ct_pers = 0
    for person in personsId :
        
        data.append([])
#        idPerson.append(person)
        
        sql = "SELECT DISTINCT traj FROM preproc WHERE person = {} AND NOT datetime = '' GROUP BY traj".format(person)
        trajectories = curs2.execute(sql).fetchall()
        
        if len(trajectories) == 0:
            raise Exception("Error: The cache did not have the requested person ID. This is most likely because the bulk cache building method was used, which is hardcoded to only load the person IDs used in the PERCOM paper.")
        
        t_ct = 0
        for t in trajectories:
            traj = t[0]
            
            sql = "SELECT idxPix, datetime FROM preproc WHERE person = {} AND traj = {} ORDER BY datetime".format(person,traj)
            rows = curs3.execute(sql)
            
            data[-1].append(list(rows))#map(lambda x: x[0],rows))
#            data[-1].append(-1)
            t_ct += 1
        
#        data[-1].pop()
#         print "person {} : {} trajectories processed.".format(person,t_ct)
        ct_pers += 1
    
    if(withDates):
        return data
    
    else:
        rtn = []
        for pers in data:
            rtn.append([])
            for traj in pers :
                rtn[-1].extend(map(lambda x: x[0],traj))
            
    print "Nb persons loaded : {}".format(ct_pers)
    print "Data loaded"
    return np.array(rtn), personsId



#Computation of Nside, corresponding to the spatial resolution :
def ComputeNside(spatialRes):
    i = 0
    listRes = []
    listNpix = []
    while True :
        listNpix.append(hp.nside2npix(2**i))
        listRes.append(510072000000000//listNpix[i])
    #     print "{} : Npix = {} ; approxRes = {}".format(i,listNpix[-1],listRes[-1]) 
        if listRes[-1] < spatialRes :
            break
        i += 1
    
    #find the nearest value :
    listRes = np.asarray(listRes)#Convertion into an numpy array :
    idx = (np.abs(listRes-spatialRes)).argmin()
    nearestSpatialRes = listRes[idx]
    npix = listNpix[idx]
    nside = hp.npix2nside(npix)
    print "nearestSpatialRes = {}, Npix = {} and Nside = {}".format(nearestSpatialRes, npix, nside)
    return nside

def buildPreprocessingTable(spatialRes,temporalRes,nest = True, personsIds = [0, 1, 2, 3, 4, 5, 7, 9, 12, 13, 14, 15, 16, 17, 22, 24, 153, 28, 30, 35, 36, 38, 39, 40, 43, 44, 50, 179, 52, 55, 68, 71, 82, 84, 85, 92, 96, 101, 104, 167, 119, 126]):

    data = []
    nside = ComputeNside(spatialRes)
   
    #connection to the DataBase :
    connectionOrig=apsw.Connection(main_geolifedb)
    
    #loading it in the memory :
    writingConn=apsw.Connection(":memory:")
#     with conn.backup("main", connectionOrig, "main") as backup:
#         backup.step() # copy whole database in one go
    
    
    curs1 = connectionOrig.cursor()
    curs2 = connectionOrig.cursor()
    curs3 = connectionOrig.cursor()
    writingCurs = writingConn.cursor()
    
    #Create table :
#     sql = "DROP TABLE IF EXISTS {}".format(tableName)
#     writingCurs.execute(sql)
    writingCurs.execute("CREATE TABLE preproc (person INT, traj INT, idxPix INT, datetime TEXT)")
    
    def getIdxPix(longitude, latitude):
        return hp.ang2pix(nside, (90 - latitude) * np.pi / 180, longitude * np.pi / 180, nest)
    
    fmt = '%Y-%m-%d %H:%M:%S'
    
    if isinstance(personsIds, str) and personsIds.lower() == 'all':
        sql = "SELECT distinct person FROM geolife  GROUP BY person "
#        WHERE person = 1
        personsIds = [ x[0] for x in curs3.execute(sql) ]
    #   
#    
    for pId in personsIds :
        
        print 'Considering s: {} t:{} Person {}'.format(spatialRes,temporalRes,pId)
        
        person = pId
#        person =0
        sql = "SELECT distinct traj FROM geolife WHERE person = {} AND NOT datetime = '' GROUP BY traj".format(person)
    #    AND traj = 20090426211055
        trajectories = curs2.execute(sql)
        
        
        prev = 0
        t_ct = 0
        #for each trajectories
        for t in trajectories:
            #select the first element of the tuple given by sqlite:
            traj = t[0]
            #make a sql query for 
            sql = "SELECT longitude ,latitude , datetime FROM geolife WHERE person = {} AND NOT datetime = '' AND traj = {} ORDER BY datetime".format(person,traj)
            rows = curs1.execute(sql)
            
            
            
            ct = 0
            points = []
            for row in rows:
                actualTime = dt.strptime(row[2],fmt)
                
                #If this is the first trajectory, set the time origin :
                if ct==0:
                    #Next time is the next time a location have to be taken
                    if t_ct == 0 :
                        #If it's the first trajectory of the person, just set the time origin at the first point
                        nextTime = actualTime
                        points.append((person,traj,int(getIdxPix(row[0],row[1])),nextTime.strftime(fmt)))
                        nextTime += temporalRes
                        
                    else :
                        #Determine the number of symbols missing. 
                        nb_loc_missing = int((actualTime- nextTime).total_seconds()//temporalRes.total_seconds()+1)
                        #Determine the next date where a location have to be taken:
                        nextTime = nb_loc_missing * temporalRes + nextTime
    #                    #If the ending and beginning locations of a gap in the GPS records are the same, 
    #                        #the user is taken as dwelling at the same location during that time.
    #                    if (int(getIdxPix(row[0],row[1])) == data[-1][-1][1]):
    #                        print "Fill blank-----------------------------------------------------", nb_loc_missing
                #Else if the time is over the next time, this mean that we have to record the point  
                elif actualTime > nextTime :
                    #Check witch point is the closer :
                    #This one: 
                    if abs(nextTime - actualTime) < abs(nextTime - dt.strptime(prev[2],fmt)):
                        #record the point :
                        points.append((person,traj,int(getIdxPix(row[0],row[1])),nextTime.strftime(fmt)))
                    #Or the previous one ?
                    else :
                        #record the point :
                        points.append((person,traj,int(getIdxPix(prev[0],prev[1])),nextTime.strftime(fmt)))
                    nextTime += temporalRes
                
                    
                ct += 1
                prev = row
            t_ct += 1
            if len(points) > 0 :
                data.append(points)            
                if len(points) > 1 :
                    writingCurs.executemany("INSERT INTO preproc (person,traj,idxPix,datetime) VALUES (?,?,?,?)",points)
                
            print "S: {} T: {} Person {}: {} trajectories processed".format( spatialRes,temporalRes, person, t_ct )

    #writing the informations about the sample :
    writingCurs.execute("CREATE TABLE infoSample (nside INT)")
    sql = "INSERT INTO infoSample VALUES ({})".format(nside)
    writingCurs.execute(sql)
    
    #Creating index :
    print "Creating index..."
    sql = """CREATE INDEX idx_person_preproc ON preproc(person);
         CREATE INDEX idx_person_traj_preproc ON preproc(person,traj);
         CREATE INDEX idx_traj_preproc ON preproc(traj);
         CREATE INDEX idx_traj_datetime_preproc ON preproc(traj,datetime);"""
    writingCurs.execute(sql)
   
   
    print "Cleaning up..."
    writingCurs.execute("vacuum")
    writingCurs.close()
   
    print "Writing out the database file..."
    
    
    #Create the database file 
    filename = "S{}T{}.sqlite".format(spatialRes, temporalRes)
    path = preprocessing_dir + "/" + filename
    ensure_dir(path)
    f = open(path, 'w')
    f.close()   
    # Now write out the database back to a file in one go
    

    connection=apsw.Connection(path)
    with connection.backup("main", writingConn, "main") as backup:
        backup.step() # copy whole database in one go

    print "Done"



def build_main_db():
    
    ensure_dir(main_geolifedb)
    
    connectionOrig=apsw.Connection(main_geolifedb)
    curs1orig = connectionOrig.cursor()
    create_table_SQL = 'CREATE TABLE geolife(person INT,  traj INT,  latitude REAL,  longitude REAL,  datetime TEXT);'
    
    curs1orig.execute( create_table_SQL )
    
    
    mainDirectory = geolife_zipfile_data_dir
    persons = os.listdir(mainDirectory)
    
    
    for person in persons:
    
        print "person n " + person
        files = os.listdir(mainDirectory + '/' + person + '/Trajectory')
        
        sqlList = []
        for trajectoryFile in files:
            
            id_trajectory = trajectoryFile.split('.')[0]
            print "\t" + id_trajectory
            with open(mainDirectory + '/' + person + '/Trajectory/' + trajectoryFile, 'r') as f:
                ct = 0
                for line in f:
                    if ct >= 6 :
                        data = line.split(',')           
                        sqlList.append((person,id_trajectory,data[0],data[1],data[-2] + " " + data[-1].replace('\n','').replace('\r','')))
                        
                    ct += 1
                          
                
            f.close()
    
        curs1orig.execute('PRAGMA journal_mode = OFF; ') # turn of journalling for speed
        curs1orig.execute('BEGIN') # this will disable autocommit for speed, bundling it into a single commit
        curs1orig.executemany('INSERT INTO geolife VALUES(?,?,?,?,?)', sqlList)  
        curs1orig.execute('END')
    
    connectionOrig.close()
    print "Done"


if __name__ == '__main__':
    def parse_timedelta(time_str):
        """
        Helper function to convert human readable time offsets to Python timedelta objects.
        This is simply to enable readability in the parameter specification.
        
        :param time_str: A time offset string in the form 'h:mm:ss'
        :type time_str: str
        """
        split = time_str.split(':')
        return timedelta(hours = int(split[0]), minutes = int(split[1]), seconds = int(split[2]))
    
    listSpatialRes = [1000,5000,10000,50000,100000,1000000,3000000,10000000,50000000]
    listRealSpatialRes = [618,2474,9896,39586,158347,633388,2533555,10134220,40536880]

    dict_SR = dict(zip(listSpatialRes,listRealSpatialRes))
    listTemporalRes = ['0:05:00','0:10:00','0:15:00','0:30:00','0:45:00','1:00:00']

    listTemporalRes = [parse_timedelta(temporalRes) for temporalRes in listTemporalRes]
    listTemporalResSecond = [temporalRes.total_seconds() for temporalRes in listTemporalRes]
    bulk_build_resolution_cache(listSpatialRes, listTemporalRes)