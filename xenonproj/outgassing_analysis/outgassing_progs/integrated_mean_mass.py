# This  program hopes to make use of the data from the RGA files to calculate an integral from the data. This should work for any given RGA data file. This program will be formatted specifically to use the .dat files in the outgassing_analysis/rga_measurements folder. This program will simply sum all of the quadrilateral components of the data to form a potential standard of analytical comparison with previous data which relies on cruder approximations.

import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ROOT as root
from scipy import integrate

file = '/home/z/xenonproj/outgassing_analysis/rga_measurements/'+'bg0714.dat' 

### The Whole introductory underlying mess is for the sake of cleaning up these awfully ugly RGA .dat files. My they really reek of 'where the hell did this come from!'. Anyways, without further ado, here we are in the end with our beautiful two arrays, our masses (in AMUs) and our pressures (in Torricellis). Magnifique.

amus = []
pressures = []

with open(file) as f:
	counter = 0
	content = f.read().splitlines()[1:]
	while counter < len(content):
		a = content[counter]
		temparray = a.split()
		del temparray[0], temparray[0]
		if len(temparray) < 3:		
			amus.append(float(temparray[0]))
			pressures.append(float(temparray[1]))
		else:
			del temparray[0]
			amus.append(int(temparray[0]))
			pressures.append(float(temparray[1]))
		counter += 1

### Here I clean up the data by removing all elements of the data with absolute value less than .18
count = 0
while count < len(pressures):
	if abs(pressures[count]) < .18e-6:
		pressures[count] = 0
	count += 1

### Now I get to take part in a solemn exercise of making out of these data points, sacrosanct rectangles, the rectangles from which Newton and Leibniz created their integral calculus, if it was in fact them who did done it. Then I get to add up all the area under these rectangles and compare them finally with those incy wincy AMUs calculated in Joey's code. Let's go find out shall we?

### This portion just finds the total area
TotalArea = 0 		### Here my friend will be our magic number for the total area under the curve
secondcounter = 0
rectanglewidth = 0.05
while secondcounter < 4000-2:
	area =  rectanglewidth* ((pressures[secondcounter+1] + pressures[secondcounter])/2)
	TotalArea = TotalArea + area
	secondcounter += 1
print(TotalArea)

### Now that the data is all cleaned up and the integral found: it's time to discriminate the peaks. I love discrimination. In so doing I will also integrate the peaks relative to each other. 



### This portion uses the total area to find the probability of hitting each amount of AMU. Using this probability for hitting on a specific AMU you can 
mean = 0.5 ### I start with 0.5 because the first 10 data points are binned in the first bin as occurs to be the case in the data file. After this first bin all the other bins are spaced apart by 20 units. All this effectively does is shift the mean value by 0.5 AMU in the positive direction.
thirdcounter = 0
while thirdcounter < 4000-2:
	area =  rectanglewidth* ((pressures[thirdcounter+1] + pressures[thirdcounter])/2)
	probability = (area/TotalArea)
	meanportion = (rectanglewidth* amus[thirdcounter]) * probability
	mean = mean + meanportion
	thirdcounter += 1
print(mean)

#### In this portion of the pretty little code we are going to do something quite magical and which lies at the heart of all scientific discourse, the lovely and illusive: graph. Plotting the points will be easy yes, and seeing their curves easier yet too. But looking deep into the curls to understand their meaning, my, that be easy pretty too.

c1 = root.TCanvas("RGA scan graph", "RGA Scan" )
g = root.TGraph()
for i in xrange(0, len(amus)):
    g.SetPoint(i, amus[i]/20.0, pressures[i])
g.SetTitle("RGA Recording for Vacuum Chamber")
g.GetXaxis().SetTitle("Mass (AMU)")
g.GetYaxis().SetTitle("Pressure (Torr)")
g.SetMarkerStyle(20)
g.SetMarkerSize(.2)
g.Draw()

### In this portion of the code, we will do something marvelous and create a smoothed out version of our data points. With this we will then be able to integrate under the graph to our heart's content.

#gs = root.TGraphSmooth()
#gs.Integral(0,200)
#grout = gs.SmoothKern(g, "normal",0.05)
#grout.Draw()
##1,"RGA Recording for Vacuum Chamber", "Pressure(Torr)", "Mass(AMU)")


### In this portion of the code we will embark on something yet too interesting: a comparison test. Here we will make use of scipy's capabilities and use built-in integration functions in the mathematics module to integrate the graph and see what it pops out. This will be used to compare both Joey's and mine method.

### First convert the above amu and pressure arrays to numpy arrays
amus_np = np.array(amus)
pressures_np = np.array(pressures)
simptergral = integrate.simps(pressures_np,amus_np,0.05)
#print(simptergral)
cumtrapzintegral = integrate.cumtrapz(pressures_np,amus_np)
#print(cumtrapzintegral[-1])

### In this portion of the coed I am going to use the given values from the Simpsons method of integration and separately with the cumulative trapezoidal means of integration find the mean mass and see how well they line up with both mine and Joey's methods for integration.
simpsmean = 0.5
fourthcounter = 0
while fourthcounter < 4000-2:
        area =  rectanglewidth* ((pressures[fourthcounter+1] + pressures[fourthcounter])/2)
        probability = (area/simptergral)
        meanportion = (rectanglewidth* amus[fourthcounter]) * probability
        simpsmean = simpsmean + meanportion
        fourthcounter += 1
#print(simpsmean)

cumtrapzmean = np.average(amus_np,0, pressures_np)/20 + 0.5
#print(amus_np)
#fifthcounter = 0
#probs = []
#while fifthcounter < 4000-2:
#	area =  rectanglewidth* ((pressures[fifthcounter+1] + pressures[fifthcounter])/2)
#	x = area/cumtrapzintegral[-1]
#	probs.append(x)
#	meanportion = (rectanglewidth* amus[fifthcounter]) * probs[fifthcounter]
#        cumtrapzmean = cumtrapzmean + meanportion
#	sumprobs = np.sum(probs)
#        fifthcounter += 1

print(cumtrapzmean)
#print( "Sum of each Independent probability is: " + str(sumprobs))

### Perhaps it's time to specify our analysis to figure out what each of these peaks really mean! In this section of the code we will integrate the main peaks and correlate them with their corresponding compounds. Then we will compare the masses of each peak to the total mass of the system.

### Check for spaces with value 0 between the peaks, when these spaces are found, the start and end of the peak is saved. Then we integrate over the peak between the two endpoints. These peak integrals are saved to a list. These peak integrals are compared to the total integral and a percentage is input in the table. These lists are then corresponded to a compound if they represent the main peak of a substance. The Substances and their percantages are printed in the output.

#integrated = []

#while fifthcounter < 4000-2:
#	if pressures[fifthcounter] == 0:
#		fifthcounter += 1
#	else:
#		while integralcounter >= fifthcounter:
#			if pressures[integralcounter] == 0:
#				break 
#			else: 
#				peakerintegration
#				fifthcounter = integralcounter + 1
#			break
#		integrated[fifthcounter]
#	return integrated 
#
