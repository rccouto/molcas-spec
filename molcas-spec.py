#!/usr/bin/env python3

import math, re, optparse, operator, os, glob
import numpy as np

def get_energy(file, level):
    E=[]
    c=0
    
    energy=0.0
    for i in open( file ).readlines():
        
        if level == "RASSCF":
            if re.search(r"::    RASSCF root number", i) is not None: # Find energy in .log
                words = i.split()
                energy = float( words[7] )  # Energy is the sixth word
                E.append(energy)
                    
        elif level == "CASPT2":
            if re.search(r"::    CASPT2 Root", i) is not None: # Find energy in .log
                words = i.split()
                energy = float( words[6] )  # Energy is the sixth word
                E.append(energy)
                    
        elif level == "MS-CASPT2":
            if re.search(r":    MS-CASPT2 Root", i) is not None: # Find energy in .log
                words = i.split()
                energy = float( words[6] )  # Energy is the sixth word
                E.append(energy)
                
        elif level == "RASSI":
            if re.search(r"::    RASSI State ", i) is not None: # Find energy in .log
                words = i.split()
                energy = float( words[6] )  # Energy is the sixth word
                E.append(energy)
        else:
            print("You forgot something, right? Maybe... level of calculation???")
            exit()
            
    if energy == 0.0:
        E.append(energy)
            
    c=c+1
    if c == 1:
        nstates=len(E)
    return E, nstates


def get_distance(atom, crd):
    C = []
    amp=[]
    
    # Get coordinates from MOLCAS files
    files=sorted(glob.iglob('*.log')) # Used to get all files in numerical order
    for file in files:
        coords = re.compile(r'\s(-?\d?\.\d+)\s+(-?\d?\.\d+)\s+(-?\d?\.\d+)\s+(\w\w\w\w\w\w\w\w) *$')
        for i in open( file ).readlines():
            if coords.search(i):
                x = float(coords.search(i).group(1) ) # Get the X Coordinate
                C.append(x)
                y = float(coords.search(i).group(2) ) # Get the Y Coordinate
                C.append(y)
                z = float(coords.search(i).group(3) ) # Get the Z Coordinate
                C.append(z)    # Save the coordinates in a List
        amp.append(C[atom*3+crd])
        del C[:]
    return amp



def get_oscillator(file):
    O=[]
    data = re.compile(r'(\d+)\s+(\d+)\s+(\d+\.\d*(?:[Ee]-?\d+)?)')
    flag=False
    with open(file, "r") as log:
        for line in log:
            if line.startswith("++ Velocity transition strengths"):
                flag=True
            if flag:
                if data.search(line):
                    if data.search(line):
                        i = float(data.search(line).group(1))   # initial state
                        O.append(i)
                        f = float(data.search(line).group(2))   # final state
                        O.append(f)
                        osc = float(data.search(line).group(3)) # oscillator strength
                        O.append(osc)
                if line.startswith("++ Length and velocity gauge comparison (spin-free states):"):
                    break
    log.close()
    states=np.array(O).reshape(-1,3)
    return states


def lorentzian(x,y,gamma,initial,final,n=1000):
    import numpy
    step=(x[-1]-x[0])/n
    xi=numpy.arange(initial,final,step)
    yi=numpy.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(y)):
            yi[i] = yi[i] + y[k] * gamma**2 / ( (xi[i]-x[k])**2 + gamma**2 )
    return xi,yi

def gaussian(x,y,sigma,initial,final,n=1000):
    import numpy
    step=(x[-1]-x[0])/n
    xi=numpy.arange(initial,final,step)
    yi=numpy.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(y)):
            yi[i] = yi[i] + y[k] * numpy.exp( -(xi[i]-x[k])**2 / sigma**2 )
    return xi,yi



# MAIN PROGRAM
def main():
    import sys
    f = optparse.OptionParser()
    # Get Type of Run
    f.add_option('-t', '--typ', type = str, default = 'help')
    # Get Number of States
    f.add_option('-s', '--st' , type = int, default = 1)    
    # Get Number of Atom from Amplitude Calculation (Only for typ=get)
    f.add_option('-a', '--atm' , type = int, default = 1)
    # Coordinate for diatomic module
    f.add_option( '-k', '--crd' , type = int, default = '0')
    # Level of calculation
    f.add_option( '-l', '--level' , type = str, default = None)
    # OpenMolcas log file
    f.add_option( '-f', '--file' , type = str, default = 'all')
    # Prinf level
    f.add_option( '-p', '--print' , type = str, default = 'list')
    # Convolution functions
    f.add_option( '-c', '--conv' , type = str, default = 'gaussian')
    # File name
    f.add_option( '-n', '--name' , type = str, default = 'spec')
    # File name
    f.add_option( '-i', '--init' , type = int, default = '1')   
    (arg, args) = f.parse_args(sys.argv[1:])

    if arg.typ == 'pec':
    
        # Get energy from .log Files
        energy,nstates=get_energy(arg.level)

        # Get bond distance
        amp = get_distance(arg.atm, arg.crd)
        
        for i in range(len(amp)):
            print (amp[i], end=' ')
   
            for j in range(nstates):
                print(energy[(nstates*i)+j], end=' ')
                print(' ')


    # Get spectrum
    elif arg.typ == 'spec':

        borh2ev=27.211385
        
        if arg.file == 'all':
            logs=sorted(glob.iglob('*.log'))
        else:
            logs=arg.file
            
        for i in range(len(logs)):
            trans=[]
            states= get_oscillator(logs[i])
            energy,nstates=get_energy(logs[i], "RASSI")

            if arg.print == "list" or arg.print == "all":
                name=arg.name
                output = "{a}-{b}-list.dat".format(a=name, b=logs[i].replace(".log", ""))
                out=open(output, "w")
                            
                out.write("# E(eV) Osc. Str. (a.u.) -- %s - List of transitions \n" % logs[i])
                for j in range(len(states)):
                    init=states[j][0].astype(int)
                    final=states[j][1].astype(int)

                    if init == arg.init:
                        out.write(" %.4f  %E \n" % ((energy[final-1]-energy[init-1])*borh2ev, states[j][2]))
                    elif arg.init == "00":
                        out.write(" %.4f  %E \n" % ((energy[final-1]-energy[init-1])*borh2ev, states[j][2]))
                out.close()
                
            if arg.print == "bars" or arg.print == "all":
                name=arg.name
                output = "{a}-{b}-bars.dat".format(a=name, b=logs[i].replace(".log", ""))
                out=open(output, "w")
                
                out.write("# E(eV) Osc. Str. (a.u.) -- %s, Bars spectrum \n" % logs[i])
                for j in range(len(states)):
                    init=states[j][0].astype(int)
                    final=states[j][1].astype(int)

                    if init == arg.init:
                        out.write(" %.4f  0.0000 \n" % ((energy[final-1]-energy[init-1])*borh2ev) )
                        out.write(" %.4f  %E \n" % ((energy[final-1]-energy[init-1])*borh2ev, states[j][2]))
                        out.write(" %.4f  0.0000 \n" % ((energy[final-1]-energy[init-1])*borh2ev) )
                out.close()
                
            if arg.print == "curve" or arg.print == "all":
                for j in range(len(states)):
                    init=states[j][0].astype(int)
                    final=states[j][1].astype(int)
                    if init == arg.init:
                        trans.append((energy[final-1]-energy[init-1])*borh2ev)
                        trans.append(states[j][2])
                trans=np.array(trans).reshape(-1,2)
                
                if arg.conv == "lorentzian":
                    name=arg.name
                    output = "{a}-{b}-curve-loren.dat".format(a=name, b=logs[i].replace(".log", ""))
                    out=open(output, "w")
                    
                    gamma = 0.1240839
                    x,y = lorentzian(trans[:,0], trans[:,1], gamma, trans[0,0]-10, trans[0,0]+10 )
                    out.write("# E(eV) Osc. Str. (a.u.) -- %s, Lorentzian Curve spectrum \n" % logs[i])
                    
                    for i in range(len(x)):
                        out.write("%s %s\n" % (x[i], y[i]))
                    out.close()
                        
                elif arg.conv == "gaussian":
                    name=arg.name
                    output = "{a}-{b}-curve-gauss.dat".format(a=name, b=logs[i].replace(".log", ""))
                    out=open(output, "w")
                        
                    sigma = 0.2
                    x,y = gaussian(trans[:,0], trans[:,1], sigma, trans[0,0]-10, trans[0,0]+10 )
                    out.write("# E(eV) Osc. Str. (a.u.) -- %s, Gaussian Curve spectrum \n" % logs[i])
                    for i in range(len(x)):
                        out.write("%s %s\n" % (x[i], y[i]))     
                    out.close()
        

        
    elif arg.typ == 'help':
        print("Do you know what do you wanna do? Apparently not. Choose the type of analysis, please!")
    
                
if __name__=="__main__":
    main()
