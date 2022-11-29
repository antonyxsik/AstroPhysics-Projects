import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits



def plotspectrum(name):
    filename = name
    data = np.genfromtxt(filename, skip_header = 17, skip_footer=1)
    pos = data[:,0]
    intensity = data[:,1]
    %matplotlib inline
    plt.figure(1)
    plt.plot(pos, intensity, 'o-', markersize = 2)
    plt.xlabel("Pixel Position")
    plt.ylabel("Intensity")
    plt.show()



def centroid_finder(x,y,bias): #Take the x and y values of the spectrum.
    s = []
    I = []
    cmean = []
    cstddev = []
    i = 0
    while i < len(x): #Goes through all of the pixels
        if y[i] > bias: #Skip over any pixel lower than the noise
            check = True
        else:
            check = False
        if check:
            s.append(x[i])
            I.append(y[i])
        if y[i-1] >= y[i] and y[i] < y[i + 1] or y[i + 1] <= bias:
#If new emission line is found, find standard deviation and mean and append it
            m = np.sum((np.multiply(s,I)) / np.sum(I))
            cmean.append(m)
#standard deviations
            std = np.sum(np.multiply(I,(s-m)**2)) / np.sum(I)
            cstddev.append(std)
#Clear the arrays
            s = []
            I = []
        i += 1 #if there are more pixels, check again.
    return [cmean,cstddev] #return the arrays of centroids and standard deviations as a new array


def linfit(x,y):
    a = np.array([[np.sum(np.power(x,2)),np.sum(x)],[np.sum(x),x.size]])
    b = np.array([[np.sum(np.multiply(x,y))],[np.sum(y)]])
    c = np.matmul(np.linalg.inv(a),b)
    return c

def linfit_error(x,y,m,c):
    N = x.size
    sig2 = 1/(N-2)*(np.sum(np.power(y-(m*x+c),2)))
    sign2=N*sig2/(N*np.sum(np.power(x,2))-np.power(np.sum(x),2))
    sigc2 = sig2*np.sum(np.power(x,2))/(N*np.sum(np.power(x,2))-np.power(np.
    ,→sum(x),2))
    return [np.sqrt(sign2),np.sqrt(sigc2)]


def reduce(science,flat,bias):
    bias = np.array(bias)/1.00001
    flat = np.array(flat)
    science = np.array(science)
    flat_norm = (flat-bias/np.median(flat-bias))
    science_reduced = science-bias/flat_norm
    return science_reduced


filenames = ["Hydrogen01870.txt", "groupc-Helium02190.txt",␣ "groupc-mercury01590.txt","groupc-neon01250.txt", "groupc-lamp01690.txt", ]
i = 0
while i<len(filenames):
    %matplotlib inline
    print("Data:", filenames[i])
    plotspectrum(filenames[i])
    i+=1



 i = 0
centroiddata = []
centroidstds = []
while i < len(filenames):
    data = np.genfromtxt(filenames[i], skip_header = 17, skip_footer = 1)
    pos = data[:,0]
    intensity = data[:,1]
    c = centroid_finder(pos, intensity, 220)
    print("file :", filenames[i])
    print(c)
    centroiddata.append(c[0])
    centroidstds.append(c[1])
    i += 1



# Wavelengths taken from the webpage recommended by professor
# Direct Link: https://physics.nist.gov/PhysRefData/ASD/lines_form.html
#Chosen Wavelengths
H_chosen_WL = np.array([410.1734, 434.0472, 486.135, 656.279])
HE_chosen_WL = np.array([257.76, 388.8648, 402.61914, 447.14802, 706.5190, 1091.,292])
HG_chosen_WL = np.array([82.1329, 167.2537, 249.2064, 567.7105, 952.0198, 1128.,71])
NE_chosen_WL = np.array([640.22480, 692.94672, 837.76070, 865.43828, 966.54200,,1117.75246])
#It is clear that the lamp has a continuous emission spectrum
#so there is no point in matching wavelengths to centroids for this particular,object
#Chosen Centroids
H_chosen_CDS = np.array([411.1832692307692, 435.0827676322418, 487.,36251279960135, 654.266424352567])
HE_chosen_CDS = np.array([265.96989528795814, 353.9975686846584, 411.,2439817826936, 445.4630620985011, 677.7912118217544,1085.5953213817227])
HG_chosen_CDS = np.array([69.48913043478261, 174.50987432675043, 258.,10568181818184, 561.0188834154351, 992.7629805786762,1218.9766777724892])
NE_chosen_CDS = np.array([645.0, 696.9963768115942, 827.6161117961366, 865.,846017699115, 981.4306358381504, 1131.9761570827488])
chosenWL = [H_chosen_WL, HE_chosen_WL, HG_chosen_WL, NE_chosen_WL]
chosenCDS = [H_chosen_CDS, HE_chosen_CDS, HG_chosen_CDS, NE_chosen_CDS]
   

names = ["Hydrogen", "Helium", "Mercury", "Neon"]

i = 0
while i<=3:
    print(names[i],"Fit")
    fit = linfit(chosenWL[i],chosenCDS[i])
    fite = linfit_error(chosenWL[i],chosenCDS[i],fit[0],fit[1])
    print("Slope, Intercept:", np.transpose(fit))
    print("Error:", fite)
    plt.figure()
    %matplotlib inline
    plt.plot(chosenWL[i], chosenCDS[i], 'o-', markersize = 2)
    x = np.array(range(int(np.min(chosenWL[i])),))
    plt.plot(x,fit[0]*x+fit[1],'k',alpha = 0.5)
    plt.xlabel("Centroid (pixel position)")
    plt.ylabel("Wavelength (nm)")
    plt.show()
    i+=1


def reduce(science,flat,bias):
    bias = np.array(bias)/1.00001
    flat = np.array(flat)
    science = np.array(science)
    flat_norm = (flat-bias/np.median(flat-bias))
    science_reduced = (science-bias)/flat_norm
    return science_reduced


Galaxy1_name = "data-2013-10-26-shane-public/b151.fits"
Galaxy2_name = "data-2013-10-26-shane-public/b158.fits"
Galaxy1 = fits.getdata(Galaxy1_name)
Galaxy2 = fits.getdata(Galaxy2_name)

bias_name = "data-2013-10-26-shane-public/b101.fits"
bias = fits.getdata(bias_name)
domeflat4_name = "data-2013-10-26-shane-public/b112.fits"
domeflat2_name = "data-2013-10-26-shane-public/b121.fits"
dome4s = fits.getdata(domeflat4_name)
dome2s = fits.getdata(domeflat2_name)
w=2020
h=20
rGalaxy1 = reduce(Galaxy1[h:279,0:w],dome4s[h:279,0:w],bias[h:279,0:w])
rGalaxy2 = reduce(Galaxy2[h:279,0:w],dome4s[h:279,0:w],bias[h:279,0:w])
%matplotlib inline


plt.figure()
plt.title("Bias Frame")
plt.imshow(bias[h:279,0:w],origin = 'lower', interpolation = 'nearest', cmap =,'pink')
plt.figure()
plt.title("Dome Frame 2s")
plt.imshow(dome2s[h:279,0:w],origin = 'lower', interpolation = 'nearest', cmap,= 'pink')
plt.figure()
plt.title("Dome Frame 4s")
plt.imshow(dome4s[h:279,0:w],origin = 'lower', interpolation = 'nearest', cmap,= 'pink')
plt.figure()
plt.title("Zw 229-015 Seyfert 1 Galaxy")
plt.imshow(Galaxy1,origin = 'lower', interpolation = 'nearest', cmap = 'pink')
plt.figure()
plt.title("Zw 229-015 Seyfert 1 Galaxy REDUCED")
plt.imshow(rGalaxy1,origin = 'lower', interpolation = 'nearest', cmap = 'pink')
plt.figure()
plt.title("3C079 Seyfert 2 Galaxy")
plt.imshow(Galaxy2,origin = 'lower', interpolation = 'nearest', cmap = 'pink')
plt.figure()
plt.title("3C079 Seyfert 2 Galaxy REDUCED")
plt.imshow(rGalaxy2,origin = 'lower', interpolation = 'nearest', cmap = 'pink')



def linearize(science,max1,max2):
    i = 0
    size = science[0].size
    spectrum = []
    m = (max2-max1)/-size
    while i<size:
        spectrum.append(science[int(max1-m*i),i])
        i+=1
    return np.array(spectrum)


%matplotlib inline
s1 = linearize(Galaxy1[h:279,0:w]-(bias[h:279,0:w]/1.000001),177,166)
plt.figure()
x = np.array(range(0,2020))
plt.plot(x, s1, 'o-', markersize = 2)
plt.xlabel("Pixel Position")
plt.ylabel("Intensity")
plt.show()
s2 = linearize(Galaxy2[h:239,0:w]-(bias[h:239,0:w]/1.000001),177,166)
plt.figure()
x = np.array(range(0,2020))
plt.plot(x, s2, 'o-', markersize = 2)
plt.xlabel("Pixel Position")
plt.ylabel("Intensity")
plt.show()


#Need to make this based on wavelength rather than pixe position, so we analyze ,the arc frame.
arc_name = "data-2013-10-26-shane-public/b100.fits"
arc = fits.getdata(arc_name)-(bias/1.0000001)
%matplotlib inline
plt.figure()
plt.imshow(arc, origin = 'lower', interpolation = 'nearest', cmap = 'magma')

larc = linearize(arc[h:279,0:w] - (bias[h:279,0:w]/1.000001),177,166)
plt.figure()
x = np.array(range(0,2020))
arc_centroids = centroid_finder(x, larc, 1000)
plt.plot(x, larc[0:2020], 'o-', markersize = 2)
plt.plot(arc_centroids, np.full(len(arc_centroids),40000), 'o', markersize = 2)
plt.xlabel("Pixel Position")
plt.ylabel("Intensity")
plt.show()

print(arc_centroids)

Ref_WL = np.array([326.1, 340.3, 361.05, 388.86, 404.66, 435.33, 447.15, 467.,82, 479.99, 508.58, 546.74])

arc_centroids = np.array(arc_centroids[0])

calibrate = linfit(arc_centroids, Ref_WL)
calibrate_error = linfit_error(arc_centroids, Ref_WL, calibrate[0],␣,calibrate[1])

cm = float(calibrate[0])
ci = float(calibrate[1])
print("Slope/Intercept:", np.transpose(calibrate))
print("Error:", calibrate_error)
%matplotlib inline
plt.figure()
plt.plot(Ref_WL, arc_centroids, 'o-', markersize = 2)
x = np.array(range(int(np.min(Ref_WL),)))
plt.plot(x, calibrate[0]*x + calibrate[1], 'k', alpha = 0.5)
plt.xlabel("Centroid (pixel position)")
plt.ylabel("Wavelength (nm)")
plt.show()


%matplotlib inline
x = np.array(range(0,2020))
plt.figure()
plt.plot(x*cm + ci, s1, 'o-', markersize = 2)
plt.plot(arc_centroids*cm + ci, np.full(len(arc_centroids),np.max(s1)), 'o',␣,markersize = 2)
plt.title("Zw 229-015 Seyfert 1 Galaxy FINAL")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.show
plt.figure()
plt.plot(x*cm + ci, s2, 'o-', markersize = 2)
plt.plot(arc_centroids*cm + ci, np.full(len(arc_centroids),np.max(s2)), 'o',␣,markersize = 2)
plt.title("3C079 Seyfert 2 Galaxy FINAL")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.show

