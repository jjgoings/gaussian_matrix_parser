# Extract matrices from Gauusian .log files and write them to text. 
# Copyright (C) 2013  Joshua J Goings
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import re
import numpy as np
import os
import sys

##########################################################
#
#       FUNCTIONS
#
##########################################################

def make_symmetric(matrix):
    return matrix+matrix.T - np.diag(np.diag(matrix))

def triangle_fill(matrix,istrt,iend,jstrt,count,element):
    for i in range(istrt,iend):
        for j in range(jstrt,i+1):
            matrix[i,j] = element[count]
            count += 1
    return matrix, count

def block_fill(matrix,istrt,nbf,jstrt,jend,count,element):
    for i in range(istrt,nbf):
        for j in range(jstrt,jend):
            matrix[i,j] = element[count]
            count += 1
    return matrix, count

def create_matrix(matrix,elements):
  ''' create lower triangular matrix from list of matrix elements
   indexed like so:
   	   [[0,0,0,...,0],
    	[1,2,0,...,0],
        [3,4,5,...,0]]
   nbf is number of basis functions
   elements is a list of matrix elements indexed like above, e.g.
      [0,1,2,3,...]
   Gaussian prints every 5 columns, so the mod5 accounts for this
  '''
  count = 0 	# count is our index
  # fill the main block, leaving remainder for triangle fill
  for i in range(0,nbf-nbf%5,5):
    matrix,count = triangle_fill(matrix,i,i+5,i,count,elements)
    matrix,count = block_fill(matrix,i+5,nbf,i,i+5,count,elements)
  # finish filling the last triangle bit
  matrix,count = triangle_fill(matrix,nbf-nbf%5,nbf,nbf-nbf%5,count,elements)
  return matrix


'''	
    Get basis functions,
       repulsion energy, set options 
'''
# Determine logfile 
if len(sys.argv) == 1:
    print "Enter name of folder containing .dat files"
    sys.exit(1)

g09file = sys.argv[1]

fileName, fileExtension = os.path.splitext(g09file)

if os.path.exists(fileName):
    sys.exit('Error: directory exists! Delete it and re-run gauss_parse')
else:
    os.makedirs(fileName)

np.set_printoptions(precision=2) 

# Extract the number of basis functions from G09 .log file
get_overlap = get_KE  = get_PE = get_ERI = False
logfile = open(g09file,'r')
for text in logfile:
    words = text.split()
    if all(x in words for x in ['Overlap']):
        get_overlap = True
    if all(x in words for x in ['Kinetic', 'Energy']):
        get_KE = True
    if all(x in words for x in ['Potential', 'Energy']):
        get_PE = True
    if all(x in words for x in ['Dumping','Two-Electron','integrals']):
        get_ERI = True
    if all(x in words for x in ['primitive','gaussians,','basis','functions,']):
        nbf = int(words[0])   # number basis functions 
    if all(x in words for x in ['nuclear','repulsion','energy','Hartrees.']):
        enuc = float(words[3]) # nuclear repulstion energy in Hartrees
    if all(x in words for x in ['alpha','beta','electrons']):
        nelec = int(words[0]) + int(words[3]) # number alpha elec + beta elec 
logfile.close()

'''	
	 Get and create overlap matrix
'''	

logfile = open(g09file,'r')
data = logfile.read()
overlap_matrix = np.zeros((nbf,nbf))
# grab all text between  "Overlap ***" and "*** Kinetic"
raw_overlap_string = re.findall(r'Overlap \*\*\*(.*?)\*\*\* Kinetic',data,re.DOTALL)
raw_overlap_string = raw_overlap_string[0].replace('D','E')
raw_overlap_elements = raw_overlap_string.split()
matrix_elements = []
for overlap_value in raw_overlap_elements:
    if 'E' in overlap_value:
        matrix_elements.append(overlap_value)
overlap = create_matrix(overlap_matrix,matrix_elements)
overlap = make_symmetric(overlap)
logfile.close()

'''	
	Get and create KE matrix
'''	

logfile = open(g09file,'r')
data = logfile.read()
KE_matrix = np.zeros((nbf,nbf))
# grab all text between  "Overlap ***" and "*** Kinetic"
raw_KE_string = re.findall(r'Kinetic Energy \*\*\*(.*?)Entering OneElI...',data,re.DOTALL)
raw_KE_string = raw_KE_string[0].replace('D','E')
raw_KE_elements = raw_KE_string.split()
matrix_elements = []
for KE_value in raw_KE_elements:
    if 'E' in KE_value:
        matrix_elements.append(KE_value)
KE = create_matrix(KE_matrix,matrix_elements)
KE = make_symmetric(KE)
logfile.close()

#print 'Kinetic Energy matrix: \n', KE

'''
	 Get and create PE matrix
'''

logfile = open(g09file,'r')
data = logfile.read()
PE_matrix = np.zeros((nbf,nbf))
# grab all text between  "Overlap ***" and "*** Kinetic"
raw_PE_string = re.findall(r'Potential Energy \*\*\*\*\*(.*?)\*\*\*\*\*\* Core Hamiltonian',data,re.DOTALL)
raw_PE_string = raw_PE_string[0].replace('D','E')
raw_PE_elements = raw_PE_string.split()
matrix_elements = []
for PE_value in raw_PE_elements:
    if 'E' in PE_value:
        matrix_elements.append(PE_value)
PE = create_matrix(PE_matrix,matrix_elements)
PE = make_symmetric(PE)
logfile.close()

#print 'Potential Energy matrix: \n', PE

'''
	 Get and create table of two electron integrals
'''

logfile = open(g09file,'r')
ERI_list = []
count = 0
for text in logfile:
    words = text.split()
    if 'I=' and 'J=' and 'K=' and 'L=' in words:
        ERI_list.append([int(words[1]),int(words[3]),int(words[5]),int(words[7]),float(words[9].replace('D','E'))])
ERI = np.array(ERI_list)
#print 'Electron repulsion integrals: \n', ERI
logfile.close()

'''       
          Write to file
'''       

if get_overlap == True:
    np.savetxt(fileName + '/overlap.dat',overlap,fmt='%.8e',delimiter = ' ')
if get_KE == True:
    np.savetxt(fileName + '/kinetic_energy.dat',KE,fmt='%.8e',delimiter = ' ')
if get_PE == True:
    np.savetxt(fileName + '/potential_energy.dat',PE,fmt='%.8e',delimiter = ' ')
if get_ERI == True:
    np.savetxt(fileName + '/two_electron_ints.dat',ERI,fmt='%d %d %d %d %.8f',delimiter = ' ')
np.savetxt(fileName + '/nuclear_repulsion.dat',np.array([enuc]),fmt='%.8f')
np.savetxt(fileName + '/number_basis_functions.dat',np.array([nbf]),fmt='%d')
np.savetxt(fileName + '/number_electrons.dat',np.array([nelec]),fmt='%d')
