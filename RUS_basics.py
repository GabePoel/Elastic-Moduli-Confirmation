#/usr/bin/python -tt

#These are the functions that should be useful for all methods of the RUS inverse problem (including the forward-solver itself).

import numpy as np
from numpy import pi,sqrt,cos,sin
import RUS_fortran_functions as funcs
from scipy.optimize import linear_sum_assignment

#solve the generalized eigenvalue problem and clean up the results
def forwardsolve(basis,mass,C,shape,dims):

    eigenvals,eigenvects = funcs.forwardsolver(basis,mass,C,shape,dims)   #create E and Gamma; solve the eigenvalue problem
    
    #note: eigenvalues are already in ascending order (with eigenvects in a corresponding order)
    #eigenvectors are the columns of a matrix; each is already normalized so that aTEa = 1

    numzeros = 0        #getting rid of 0Hz freqs, and corresponding eigenvects
    for ii in range(0,len(eigenvals)):
        if eigenvals[ii] < 10**(-4):   #(cutting all frequencies below a few kHz; note that eigenvals = w^2 in MHz^2)
            numzeros +=1
    eigenvals = eigenvals[numzeros:]
    eigenvects = eigenvects[:,numzeros:]

    omegas = sqrt(eigenvals).real   #(units: MHz)
    freqs = np.multiply(1000/(2*np.pi),omegas)  #(new units: kHz)

    degindex = []       #getting rid of degenerate freqs, and corresponding eigenvects
    for ii in range(0,len(freqs)-1):
        deg = close_floats(freqs[ii],freqs[ii+1],1e-2)  #(anything within 0.01 kHz counts as degenerate)
        if deg:
            degindex.append(ii)
    for jj in range(len(degindex)-1,-1,-1):   #(we have to delete working back from the end, to avoid index confusion)
        freqs = np.delete(freqs,degindex[jj])
        eigenvects = np.delete(eigenvects,degindex[jj],1)

    return freqs,eigenvects




#create a list of the functions x^l*y^m*z^n such that l+m+n <= N
def basisfuncs(N):
    numcombos = (N+1)*(N+2)*(N+3)//6
    basis = np.zeros((numcombos,3),dtype=int)
    index = 0
    for ll in range (0,N+1):
        for mm in range (0,N+1):
            for nn in range (0,N+1):
                if ll+mm+nn < N+1:
                    basis[index] = [ll,mm,nn]
                    index += 1
    return basis



def index_compression(i,j):
    if i == j:
        return i
    elif (i+j) == 3:
        return 3
    elif (i+j) == 2:
        return 4
    else:
        return 5

#a function to populate the full 3x3x3x3 tensor Cijkl from the independent Voigt-notation elements C11, C12, etc.
def build_C(indep_elems):
    
    #build C_Voigt
    #-------------

    num_indep_elems = len(indep_elems)

    #all symmetries from isotropic up to orthorhombic have nine non-zero entries in common.
    #I call these C9:
    #C9 = [c11, c22, c33, c23, c13, c12, c44, c55, c66]

    #in addition, trigonal crystals and some types of tetragonal crystal have other non-zero entries.
    #there are 8 that can be non-zero in all. by default, they are zero; the entries only become non-zero in certain cases.
    c14,c15,c16,c24,c25,c26,c46,c56 = 0,0,0,0,0,0,0,0
    

    if num_indep_elems < 4:
        #2: isotropic
        #3: cubic
        
        c11 = indep_elems[0]
        c44 = indep_elems[1]
        if num_indep_elems == 2:
            c12 = c11 - 2*c44
        else:
            c12 = indep_elems[2]
        C9 = [c11,c11,c11,c12,c12,c12,c44,c44,c44]

    if num_indep_elems > 4 and num_indep_elems < 9:
        #5: hexagonal
        #7/8: (really 6/7; final element distinguishes trigonal/tetragonal)

        
        c33 = indep_elems[0]
        c23 = indep_elems[1]
        c12 = indep_elems[2]
        c44 = indep_elems[3]
        c66 = indep_elems[4]
        if num_indep_elems == 5: #hexagonal
            c11 = c12 + 2*c66
        else:
            if indep_elems[-1] == True: #trigonal
                c11 = c12 + 2*c66
                c14 = indep_elems[5]
                c24 = -c14
                c56 = c14
                if num_indep_elems == 8:
                    c15 = indep_elems[6]
                    c25 = -c15
                    c46 = -c15
            
            else: #tetragonal
                c11 = indep_elems[5]
                if num_indep_elems == 8:
                    c16 = indep_elems[6]
                    c26 = -c16   
            
        C9 = [c11,c11,c33,c23,c23,c12,c44,c44,c66]



    if num_indep_elems == 9:
        #orthorhombic
        
        C9 = indep_elems
        
    #turn C9 and other elems into 6x6 array
    #note that this is not the FULL Voigt-notation matrix, as the lower-left triangle is all 0s
    C_Voigt =  np.array([[C9[0],C9[5],C9[4],c14,c15,c16],
                         [0,C9[1],C9[3],c24,c25,c26],
                         [0,0,C9[2],0,0,0],
                         [0,0,0,C9[6],0,c46],
                         [0,0,0,0,C9[7],c56],
                         [0,0,0,0,0,C9[8]]])

    #turning this into the full C_Voigt
    for i in range(0,6):
        for j in range(0,6):
            if i>j:
                C_Voigt[i,j] = C_Voigt[j,i]
                
    C = Voigttotensor(C_Voigt)
    return C


#take a rank-4 tensor and turn it into the true (symmetric) C_Voigt
def tensortoVoigt(C_tens):

    C_Voigt = np.zeros((6,6),dtype=np.float64)
    
    for i in [0,1,2]:
        for j in [0,1,2]:
            for k in [0,1,2]:
                for l in [0,1,2]:
                    index_one = index_compression(i,j)
                    index_two = index_compression(k,l)
                    C_Voigt[index_one,index_two]=C_tens[i,j,k,l]

    return C_Voigt

#build rank-4 tensor from C_Voigt
def Voigttotensor(C_Voigt):

    C = np.zeros((3,3,3,3),dtype=np.float64)

    for i in [0,1,2]:
        for j in [0,1,2]:
            for k in [0,1,2]:
                for l in [0,1,2]:
                    index_one = index_compression(i,j)
                    index_two = index_compression(k,l)
                    C[i,j,k,l] = C_Voigt[index_one,index_two]
    return C


#tests if two floating-point values are approximately equal
def close_floats(a,b,abs_tol):
    min = a-abs_tol
    max = a+abs_tol
    if min <= b and b <= max:
        return True #they are approximately equal
    else:
        return False #they are not approximately equal



#creates an "assignment" between the measured and calculated frequences, minimizing total distance between matched-up freqs and making sure the pairings stay "in order"
def assign(calc,meas):

    #create the "cost matrix" to use the algorithm on
    
    #this assumes len(calc) >= len(meas), always
    distances = np.zeros((len(calc),len(calc)))
    for i in range (0,len(meas)):
        for j in range (0,len(calc)):
            distances[i,j] = np.abs(meas[i]-calc[j])
            if i > j:
                distances[i,j] += calc[-1]*10000  #add a severe penalty for assigning a measurement to a lower-ordered calc
    #note: we have len(calc)-len(meas) "dummy" rows at the end of our matrix, filled with 0s, just so it's square

    #use the "Hungarian algorithm"
    row_ind, col_ind = linear_sum_assignment(distances)              
    return col_ind


#rotate counterclockwise.  CAUTION: rotations do not commute. this only gives ONE ORDER of rotations (about Z, then Y, then X)
#input C is the full 3x3x3x3 tensor
def rotate_C_general(C,thetaX,thetaY,thetaZ):
    cX = cos(thetaX)
    sX = sin(thetaX)
    cY = cos(thetaY)
    sY = sin(thetaY)
    cZ = cos(thetaZ)
    sZ = sin(thetaZ)
    RotX = np.array([[1,0,0,0,0,0],
                     [0,cX**2,sX**2,2*cX*sX,0,0],
                     [0,sX**2,cX**2,-2*cX*sX,0,0],
                     [0,-cX*sX,cX*sX,cX**2-sX**2,0,0],
                     [0,0,0,0,cX,-sX],
                     [0,0,0,0,sX,cX]])
    RotY = np.array([[cY**2,0,sY**2,0,2*cY*sY,0],
                     [0,1,0,0,0,0],
                     [sY**2,0,cY**2,0,-2*cY*sY,0],
                     [0,0,0,cY,0,-sY],
                     [-cY*sY,0,cY*sY,0,cY**2-sY**2,0],
                     [0,0,0,sY,0,cY]])
    RotZ = np.array([[cZ**2,sZ**2,0,0,0,2*cZ*sZ],
                     [sZ**2,cZ**2,0,0,0,-2*cZ*sZ],
                     [0,0,1,0,0,0],
                     [0,0,0,cZ,sZ,0],
                     [0,0,0,-sZ,cZ,0],
                     [-cZ*sZ,cZ*sZ,0,0,0,cZ**2-sZ**2]])
    RotXT = RotX.T
    RotYT = RotY.T
    RotZT = RotZ.T

    C_Voigt = tensortoVoigt(C)
    
    C_rot = RotX.dot(RotY).dot(RotZ).dot(C_Voigt).dot(RotZT).dot(RotYT).dot(RotXT)

    C_rot_tensor = Voigttotensor(C_rot)
    
    return C_rot_tensor
