'''
Created on 26 Jul 2013

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


gpuID = 0

def setGPU( id_in ):
    global gpuID
    gpuID = id_in
    
def getGPU():
    global gpuID
    return gpuID