import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

"""
In order to find the desired parameters of the etching experiment,
this code allows the user to define what is
desired and say what is necessary to be input.
"""


starting_blurb = """
_______________________________________________________________________________________________________________
Microlens Fabrication program - MiFabPro
Created by Helen Witkowska

This code is designed to give the parameters desired by the Strathclyde O2/Ar ICP etching to create microlenses
The inputs required for each calculation is different, so the calculation desired should be input:
When prompted, input the option required.

    1) Input the height of the microlens, h, the radius of the lens, r, and output the Radius of curvature of the lens
    2) Input the RoC, lens diameter to give the height of the microlens
    3) Input the RoC, height and radius of the microlens to find the thickness of photoresist required.
    4) Input the photoresist height and reflow temperature and get the maximum possible diameter of the microlens as the output
    5) Input the photoresist thickness, diameter to find the height/shape of the reflowed photoresist and the ROC, height and the minimum etching time to get the lens
    6) Input 2/3 of ROC, height and diameter of the lens, find the focal length and NA of the diamond microlens
    7) Input the desired ROC and diameter of the lens, get the minimum necessary h0, the minimum etch time, the minimum reflow temp, and the shape of the reflowed PR and lens
    8) Input the desired focal length and reflow temperature and get the required lens radius and height.
    9) Input h0, r_t, r_b, r_reflow, h_reflow, r_lens, h to get shape on graph and evaporation percentage
    10) Input the d and h0 to find the minimum reflow temperature
   
NB. There are several assumptions in place. The first being that the lens forms a perfect circle on the diamond 
      surface. The second assumption is that all spheres are perfectly smooth. The third is that every length
      measurement is in microns.
______________________________________________________________________________________________________________
"""
print(starting_blurb)
X = input('Choose an option = ')


def radius_curvature(r,h):
    #ROC calculation using ROC = (r^2 + h^2)/2h
    roc = (r**2 + h**2)/(2*h)
    return roc

def height(r,ROC):
    #Rearranging equation (y+ROC)^2 + x^2 = ROC^2 to get y = ROC +- sqrt(ROC^2 - r^2)
    #Assumes a spherical approximation of the material
    y = ROC**2
    s = r**2
    root= y- s
    #h is the difference between y(x=0) and y(x=radius_lens) so use minus bit of the +-
    if ROC>=r:   #ROC should usually be >=r otherwise need an ellipsoidal approximation.
        h= ROC - np.sqrt(root)
    else:
        print('Error, no real values for h found - need ellipsoidal approximation')
    return h

 

def microlens(ROC,r,x,h,thickness_diamond):
    #Create the outline of the microlens using spherical approximation
    y=np.zeros(len(x))
    i = 0
    for i in range(0,len(x)):
        #creates circular shape
        root = ROC**2 - (x[i])**2
        y[i] = np.sqrt(root)
        i+=1

    l = ROC - h
    j=0
    z=np.zeros(len(x))
    for j in range(0,len(x)):
        #Links together the circular segment and the flat sections of the sample
        if x[j] < -r:
            z[j] = thickness_diamond
            j+=1
        elif x[j] > r:
            z[j] = thickness_diamond
            j+=1
        else:
            z[j] = y[j] - l + thickness_diamond
            j+=1

    return z

def photoresist_profile(x,r,thickness,h):
    #This assumes a non spherical photoresist, creates an ellipsoidal profile
    #using equation x^2/r^2 +y^2/h^2 = 1
    y = np.zeros(len(x))
    z = np.zeros(len(x))
    h_sq = h**2
    r_sq = r**2
    root = np.zeros(len(x))
    for i in range(0,len(x)):
        root[i] = h_sq*(1 - (x[i]**2)/r_sq)
        z[i] = np.sqrt(root[i])

    k = 0    
    for j in range(0,len(x)):
        if x[j]<-r:
            y[j] = thickness
        elif x[j]>r:
            y[j] = thickness
        else: 
            y[j] = z[j] +thickness
            

    plt.plot(x,y)
    plt.show()
    return y
    
def growth_start(t,x,rate_resist,rate_diamond,r,y,dt,angle,dx):
    #Reverse of the etching, start with the microlens and end with the PR on the substrate
    #Radius of photoresist and diamond must be the same
    time= 0
    t_max = t/dt
    m = np.zeros(len(y))
    m[:]+= y[:]
    for time in range(0,int(t_max)):
        i=0
        #This part works out the height, m, at which the diamond started 
        for i in range(0,len(x)):
            if x[i]<-r:
                m[i]+= dt * rate_diamond
            elif x[i]>r:
                m[i]+= dt * rate_diamond
         
            i+=1    
        time+=1

    j=0
    #a and b are the time parameters for photoresist and diamond etching respectively
    
    a = np.zeros(len(x))
    b = np.zeros(len(x))
    for j in range(0,len(x)):
        if x[j]<=-r or x[j]>=r:
            #Outside the radius of the lens, all etching is diamond etching
            b[j] = t
            a[j] = 0
        else:
            #Inside the radius of the lens, there is a mix of diamond and PR etching
            #Since the diamond started at a uniform thickness, the diamond etching is
            #limited to the distance between the etch end point and the height at which the diamond ends.
            #Then the Photoresist etching must be for the rest of the time.
            b[j] = (m[0] - m[j])/rate_diamond
            a[j] = t - b[j]
   
    k= 0
    z=np.zeros(len(x))
    for k in range(0,len(x)):
        #Combine this information with the etch rates and end etch position to get the end 
        z[k] = a[k]*rate_resist + b[k]*rate_diamond + y[k]
        k+=1
        
    plt.plot(x,z)
    plt.show()
    
    total = 0
    PR = np.zeros(len(x))
    PR[:] = z[:] - z[0]

    middle = int((len(x)-1) /2)

    print ('PR at centre = ',PR[middle])

    #Integrate the profile of the photoresist to get the crossectional area of the Photoresist
    total = sum(PR[:])*dx
    print('PR ',PR)

    
    if angle == 0:
        height = total/(2*r)
    else:
        #If the exposure/development stage gives off vertical walls
        height = total/(2*r- np.sin(angle))

    photo = np.zeros(len(x))
    for i in range(0,len(x)):
        if x[i]<-r:
            photo[i] = z[0]
        elif x[i]>r:
            photo[i] = z[0]
        else:
            photo[i] = height
    return height,PR,photo


def h_height(trunc,dx,x,r,vol_0,con,h_0):
    #Input the initial height of the PR, when it was spin coated onto the substrate
    #And find the shape of the reflowed photoresist
    i=0
    #trunc is the point at which the sequence truncates and returns an answer if it has not been reached already

    vol_est = np.zeros(trunc)
    h = np.zeros(trunc +1)
    h[0] = h_0
    k=0
    for i in range(0,trunc):
        r_sq = r**2
        twoh = 2* h[i]

        ROC = (r**2)/(2*h[i]) + h[i]/2
        l = ROC - h[i]

        if ROC < r:
            h[i] = 2*h_0
            i+=1
  
        angle = (2* h[i]*r)/(r**2 + h[i]**2) 
        arc = np.arcsin(angle)
        brac1= (r**2)/2 + (r**4)/(4*(h[i]**2)) + (h[i]**2)/4
        #The estimated volume is found by geometry, the area of the sector minus the area of the triangle
        vol_est[i] = arc*brac1 - r*l
        vol_diff = vol_0 - vol_est[i]
        print('Vol difference = ',vol_diff)

        #Numerical methods of convergence x(n+1) = x(n) +- x(n)*(actual volume - estimated volume)/2*actual volume
        # The +- depends on if the differnce is positive or negative
        if i == trunc -1:
            p= h[i]
            
        if vol_diff > con :
            h[i+1] = h[i] *(1 + (vol_diff/(2*vol_0)))
            i+=1
        elif vol_diff < -con:
            h[i+1] = h[i]*(1 - (vol_diff/(2*vol_0)))
            i+=1
        #If the differnece is less than the desired convergence level, the answer has been found.
        elif vol_diff < con:
            if vol_diff > -con:
                k=i
                print('h found at iteration ',i)
                print(' h = ',h)
                p = h[i]
                i=trunc
                break
    return h ,p    


def developed(r_dev, x_grid,h0):
    
    y = np.zeros(x_range)
    for i in range(0,x_range):
        if x[i]< -r_b or x[i]>r_b:
            y[i] = 0
        elif x[i]<r_t and x[i]>-r_t:
            y[i] = h
            
    return y        
            
        



def main(X):
    #Contains all the program inputs and collections of useful functions
    if X == '1':
        #Input the height of the microlens, h, the radius of the lens, r, and output the Radius of curvature of the lens
        print('Diameter of the microlens: ')
        d_lens = int(input())
        r_lens = d_lens/2
        r_diamond = r_lens + 2

        h = float(input('Height of the microlens = '))
                              
        ROC = radius_curvature(r_lens, h)
        print('RoC = ',ROC)                      
                              
    elif X == '2':
        #Input the RoC, lens diameter to give the height of the microlens
        print('Diameter of the microlens: ')
        d_lens = int(input())
        r_lens = d_lens/2
        r_diamond = r_lens + 2

        ROC_lens = float(input('ROC of the microlens = '))
        h = height(r_lens,ROC_lens)
        print('Height of microlens = ',h)

    elif X == '3':
        #Input the RoC, height and radius of the microlens to find the thickness of photoresist required.

        G = input('Change from the standard parameters? yes/no ')
        if G == 'yes':
            
        
            I = input('Change from standard dx value (0.1)? yes/no ')
            if I == 'no':
                dx = 0.1
            else:
                dx =float( input('dx value (ensuring d_lens/dx is an integer) = '))


            J = input('Change from standard final diamond thickness (0.05)? yes/no')
            if J =='no':
                thickness = 0.05
            else:
                thickness = float(input('Thickness of diamond = '))

            #Oxygen/Argon etch rate
            DR = input('Change the diamond etch rate from 0.22 microns/min? yes/no')
            if DR == 'no':
                rate_diamond = 0.22
            else:
                rate_diamond = float(input('Diamond etch rate = '))

            PR =  input('Change the PR to diamond etch rate ratio from 0.14? yes/no')  
            if PR == 'no':
                rate_resist = rate_diamond/0.14
            else:
                rate_resist = float(input('Photoresist to Diamond etch rate ratio = '))

            ANG = input('Are the walls of the photoresist cyliders vertical? yes/no ')
            if ANG == 'yes':
                angle = 0
            else:
                angle = float(input('Angle off vertical (positive /-\, negative \-/) (radians) = '))    

        else:
            #Standard parameters
            dx = 0.1
            thickness = 0.05
            rate_diamond = 0.22
            rate_resist = rate_diamond /0.14
            angle = 0

        #Inputs 
        d_lens = float(input('Diameter of the microlens: (ensure diameter/dx = integer) '))
        r_lens = d_lens/2
        r_diamond = r_lens + 2

        ROC_lens = float(input('ROC of the microlens = '))
        h = height(r_lens,ROC_lens)
        
        x_grid_size = int((2*r_diamond/dx)+1)
        x_grid = np.linspace(-r_diamond, r_diamond,x_grid_size)
        

        t = int(input('Time taken for the etch = '))
        dt = 0.1

             
        y = microlens(ROC_lens,r_lens,x_grid,h,thickness)
        print('y = ',y)
        g,PR,photoresist = growth_start(t,x_grid,rate_resist,rate_diamond,r_lens,y,dt,angle,dx)
        
        print('Origninal thickness of photoresist = ', g)



        
        plt.plot(x_grid,PR)
        plt.plot(x_grid,y)
        plt.plot(x_grid,photoresist)
        
        I = input('FInd the percentage of PR evaporated? yes/no ')
        if I =='yes':
            J = input('Did the diameter of the PR change during reflow? yes/no')
            if J == 'yes':
                r_pillar = float(input('Original radius of PR pillar = '))
            else:
                r_pillar = r_lens
            #equating vol of sphere cap and PR pillar using : SC*X = PR
            middle = int((len(x_grid)-1) /2)
            PR_squared = PR[middle]**2
            r_squared_l = r_lens**2
            r_squared_p = r_pillar**2
            check = PR[middle] *(3*r_squared_l +PR_squared)
            u = 6*r_squared_p *g
            evap = (1- check/u)*100
            
            print('Evaporation of solvent = ',evap)

    elif X == '4':
        #Input the photoresist height and reflow temperature and get the maximum possible diameter
        #of the microlens as the output
 
        T_p = 0.0266
        T_g = 352.2

        Temp = int(input('Reflow temp = '))
        T = Temp + 273.15
        h_0 = float(input('Photoresist thickness = '))
 
        d = h_0 * ((T-T_g)/T_p)**(1/3)
        D = d*2
        print('Maximum lens diameter ',D)

    elif X == '5':
        # Input the photoresist thickness, diameter to find the height/shape of the reflowed photoresist
        # and the ROC, height and the minimum etching time to get the lens


        G = input('Change from the standard parameters? yes/no ')
        if G == 'yes':
            
        
            I = input('Change from standard dx value (0.1)? yes/no ')
            if I == 'no':
                dx = 0.1
            else:
                dx =float( input('dx value (ensuring d_lens/dx is an integer) = '))


            J = input('Change from standard final diamond thickness (0.05)? yes/no')
            if J =='no':
                thickness = 0.05
            else:
                thickness = float(input('Thickness of diamond = '))

            DR = input('Change the diamond etch rate from 0.22 microns/min? yes/no')
            if DR == 'no':
                rate_diamond = 0.22
            else:
                rate_diamond = float(input('Diamond etch rate = '))

            PR =  input('Change the PR to diamond etch rate ratio from 0.14? yes/no')  
            if PR == 'no':
                rate_resist = rate_diamond/0.14
            else:
                rate_resist = float(input('Photoresist to Diamond etch rate ratio = '))

            ANG = input('Are the walls of the photoresist cyliders vertical? yes/no ')
            if ANG == 'yes':
                angle = 0
            else:
                angle = float(input('Angle off vertical (positive /-\, negative \-/) (radians) = '))

            CON = input('Change the level of convergence of volumes from 1? yes/no')
            if CON == 'no':
                convergence = 0.001
            else:
                convergence = float(input('Convergence = '))

            TRUNC = input('Change the number of iterations before truncation from 15? yes/no')
            if TRUNC == 'no':
                trunc = 15
            else:
                trunc = float(input('Convergence = '))    

        else:
            dx = 0.1
            thickness = 0.05
            rate_diamond = 0.22
            rate_resist = rate_diamond /0.14
            angle = 0
            convergence = 0.001  
            trunc = 15         #Much greater than usually needed
            
        P = input('Input the resist cylinder diameter or lens diameter? resist/lens')
        if P =='resist':
            d_resist = float(input('Diameter of the resist cylinder: (ensure diameter/dx = integer) '))
            d_lens = 0.8 * d_resist
            print('d_lens = ',d_lens)
        else:
            d_lens = float(input('Diameter of the lens : (ensure diameter/dx integer) '))
            d_resist = d_lens/0.8
            print('d_resist = ',d_resist)
        
        r_lens = d_lens/2
        r_diamond = r_lens + 2

        E = input('Is there solvent evaporation? yes/no ')
        if E == 'yes':
            evap = float(input('Percentage of evaporation = '))
        else:
            evap = 0
                         
        
        x_grid_size = int((2*r_diamond/dx)+1)
        x_grid = np.linspace(-r_diamond, r_diamond,x_grid_size)
        

        h_0 = float(input('Photoresist thickness - ') ) 
        vol = h_0 * d_resist

        h,h_final=h_height(trunc,dx,x_grid,r_lens,vol,convergence,h_0)
        print('h values ',h)
        print('Final h ',h_final)
        print('Actual volume = ', vol)

        L = input('Find the ROC and height of the lens? ')
        if L == 'yes':
            t_pr = h_final/rate_resist     #time to etch all of the photoresist
            h_lens = rate_diamond * t_pr
            print('height of microlens = ',h_lens)
            roc =  radius_curvature(r_lens,h_lens)
            print('ROC = ',roc)

    elif X=='6':
        #Input 2/3 of ROC, height and diameter of the lens, find the focal length and NA of the diamond microlens

        I = input('Which information will you be inputing? ROC/h/d Include a space between the names and write in the order shown - ')
        if I =='ROC h':
            ROC = float(input('Input the ROC: '))
            h = float(input('Input the h: '))
            r_sq = ROC * 2 * h - h**2
            r = np.sqrt(r_sq)
            focal_length = ROC/1.42
            NA = r /focal_length
            print('focal length = ',focal_length)
            print('NA = ',NA)
            
        elif I == 'ROC d':
            ROC = float(input('Input the ROC: '))
            d = float(input('Input the diameter: '))
            r = d/2
            focal_length = ROC / 1.42
            NA = r/focal_length
            print('focal length = ',focal_length)
            print('NA = ',NA)
            
        elif I == 'h d':
            h = float(input('Input the height of the microlens: '))
            d = float(input('Input the diameter of the microlens: '))
            r = d/2
            ROC = h/2 + (r**2)/(2*h)
            focal_length = ROC/1.42
            NA = r/focal_length
            print('focal length = ',focal_length)
            print('NA = ',NA)                

    elif X=='7':
        #Input the desired ROC and diameter of the lens, get the minimum necessary h0, the minimum etch time, the minimum reflow temp, and the shape of the reflowed PR and lens

        G = input('Change from standard parameters? yes/no ')
        if G =='yes':
            
            I = input('Change from standard dx value (0.1)? yes/no ')
            if I == 'no':
                dx = 0.1
            else:
                dx =float( input('dx value (ensuring d_lens/dx is an integer) = '))

            DR = input('Change the diamond etch rate from 0.22 microns/min? yes/no')
            if DR == 'no':
                rate_diamond = 0.22
            else:
                rate_diamond = float(input('Diamond etch rate = '))

            PR =  input('Change the PR to diamond etch rate ratio from 0.14? yes/no')  
            if PR == 'no':
                rate_resist = rate_diamond/0.14
            else:
                rate_resist = float(input('Photoresist to Diamond etch rate ratio = '))

            J = input('Change from standard final diamond thickness (0.05)? yes/no')
            if J =='no':
                thickness = 0.05
            else:
                thickness = float(input('Thickness of diamond = '))

            ANG = input('Are the walls of the photoresist cyliders vertical? yes/no ')
            if ANG == 'yes':
                angle = 0
            else:
                angle = float(input('Angle off vertical (positive /-\, negative \-/) (radians) = '))      

        else:
            dx = 0.1
            rate_diamond = 0.22
            rate_resist = rate_diamond/0.14
            thickness = 0.05
            dt = 0.1
            angle = 0
            
 
             #Adjusted - need to remove this adjustion   
        #Input the desired ROC and diameter of the lens, get the minimum necessary h0 and the minimum etch time
        ROC_lens = float(input('Input the ROC of the lens= '))
        ROC_PR = float(input('Input the ROC of the photoresist = '))
        d = float(input('Input the diameter of the lens = '))
        d_PR = float(input('Input the diameter of the Photoresist = '))
        r_lens = d/2
        r_PR = d_PR/2
        r_diamond = r_lens+2
        root_lens = ROC_lens**2 - r_lens**2
        root_PR = ROC_PR**2 - r_PR**2
        h_lens = ROC_lens - np.sqrt(root_lens)
        h_PR = ROC_PR - np.sqrt(root_PR)
        t = h_lens * rate_resist
        print('Minimum etch time = ',t)
        x_grid_size = int((2*r_diamond/dx)+1)
        x_grid = np.linspace(-r_diamond, r_diamond,x_grid_size)
        y = microlens(ROC_PR,r_lens,x_grid,h,thickness)
        g,PR = growth_start(t,x_grid,rate_resist,rate_diamond,r_PR,y,dt,angle,dx)
        print('Initial photoresist thickness = ',g)
        plt.plot(x_grid,PR)
        plt.plot(x_grid,y)
        plt.show()

        cube = r_lens /g
        Temp = cube**3 * 26.6/1000 +352.2 -273
        print('Minimum reflow temperature required (celsius) = ',Temp)

        

    elif X == '8':

        #Input the desired focal length and reflow temperature to get the required lens radius and height
        I = input('Is the sample diamond? yes/no ')
        if I == 'yes':
            n = 2.42
        else:
            n = float(input('What is the reflective index of the sample? '))
            
        f = float(input('Desired focal length = '))
        '''
        J = input('Input h or r? ')
        if J == 'h':
            h = float(input('h = '))
            r = np.sqrt(h*(2.84*f-h))
        elif J == 'r':
            r = float(input('r = '))
            h = 
        '''            
        ROC = f*(n-1)
        two_roc = 2*ROC
        four = int(4*two_roc)
        h = np.linspace(0,two_roc,four)
        r = np.zeros(len(h))
        for i in range(0,len(h)):
            r[i] = np.sqrt(h[i]*((n-1)*f-h[i]))

        plt.plot(h,r)
        plt.ylabel('Radius of lens')
        plt.xlabel('height of lens')
        plt.show()


    elif X == '9':
        
        #Input h0, r_t, r_b, r_reflow, h_reflow, r_lens, h to get shape on graph and evaporation percentage

        
        G = input('Change from the standard parameters? yes/no ')
        if G == 'yes':
            
        
            I = input('Change from standard dx value (0.1)? yes/no ')
            if I == 'no':
                dx = 0.1
            else:
                dx =float( input('dx value (ensuring d_lens/dx is an integer) = '))


            J = input('Change from standard final diamond thickness (0.05)? yes/no')
            if J =='no':
                thickness = 0.05
            else:
                thickness = float(input('Thickness of diamond = '))

            DR = input('Change the diamond etch rate from 0.22 microns/min? yes/no')
            if DR == 'no':
                rate_diamond = 0.22
            else:
                rate_diamond = float(input('Diamond etch rate = '))

            PR =  input('Change the PR to diamond etch rate ratio from 0.14? yes/no')  
            if PR == 'no':
                rate_resist = rate_diamond/0.14
            else:
                rate_resist = float(input('Photoresist to Diamond etch rate ratio = '))

            CON = input('Change the level of convergence of volumes from 1? yes/no')
            if CON == 'no':
                convergence = 0.001
            else:
                convergence = float(input('Convergence = '))

            TRUNC = input('Change the number of iterations before truncation from 15? yes/no')
            if TRUNC == 'no':
                trunc = 15
            else:
                trunc = float(input('Convergence = '))    

        else:
            dx = 0.1
            thickness = 0.05
            rate_diamond = 0.22
            rate_resist = rate_diamond /0.14

            convergence = 0.001  
            trunc = 15         #Much greater than usually needed
        r_dev = float(input('Developed pillar radius = '))
        h0 = float(input('Initial thickness of photoresist = '))
        r_reflow = float(input('Radius of reflowed spherical cap = '))
        h_reflow = float(input('Height of reflowed spherical cap = '))
        r_lens = float(input('Radius of lens = '))
        h = float(input('Height of lens = '))

        radius = int(r_reflow)
        n = int((2*radius +5)/dx)
        x_grid =  np.linspace(-(radius+2),radius+2,n)

        ROC_lens = h/2 + r_lens**2/(2*h)
        y_dev = developed(r_dev, x_grid,h0)
        y_PR = photoresist_profile(x_grid,r_reflow,thickness,h_reflow)          
        y_diamond =  microlens(ROC_lens,r_lens,x_grid,h,thickness)

        plt.plot(x_grid,y_PR, label = 'Reflowed photoresist')
        plt.plot(x_grid,y_diamond, label = 'Diamond SIL')
        plt.plot(x_grid,y_dev, label = 'Developed photoresist')
        plt.legend(loc='upper right')
        plt.ylabel('Height (microns)')
        plt.xlabel('Width (microns)')
        plt.show()

    elif X == '10':
        #Input the d and h0 to find the minimum reflow temperature
        d = float(input('Radius of reflowed sphere = '))
        h0 = float(input('Initial height of photoresist = '))
        Tp = 26.6/1000
        Tg = 352.2
        T = Tp *(d/h0)**3 +Tg
        Temp = T -273.15
        print('Minimum reflow temperature = ', Temp)
        
            
    else:
        print('Error, option makes no sense, please pick from one of the options')

    
m = main(X)
        

