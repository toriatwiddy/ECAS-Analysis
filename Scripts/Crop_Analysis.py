from skimage import io, measure
from PIL import Image
import pandas as pd
import numpy as np
#from numpy import asarray
from scipy import spatial
import matplotlib.pyplot as plt

# Location to save output:

root_directory=r"D:/ECAS_Cohort/FL_Crops/Moderate"

pathList=[]

pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/1")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/2")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/3")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/4")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/5")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/6")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/7")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/8")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/9")
pathList.append(r"D:/ECAS_Cohort/FL_Crops/Moderate/SD045_15_crops/10")

pixel_size=325

# Function to open images:

def open_img(toload):

    img=io.imread(toload)

    return img

# Function to threshold images using fixed threshold value:

def threshold_img_fixed(input_img,threshold_value):

    binary_img=input_img>threshold_value

    return binary_img

# Function to label the features in the thresholded image:

def label_image(input_img):

    labelled_img=measure.label(input_img)
    features_number=labelled_img.max()

    return features_number,labelled_img

# Function to take a labelled image + the original image and measure properties:

def analyse_labelled_img(labelled_img,intensity_img):

    measure_img=measure.regionprops_table(labelled_img,intensity_img,properties=('area',
                                                                    'perimeter',
                                                                    'centroid',
                                                                    'orientation',
                                                                    'major_axis_length',
                                                                    'minor_axis_length',
                                                                    'mean_intensity',
                                                                    'max_intensity'))
    measure_df=pd.DataFrame(measure_img)

    return measure_df

# Function to analyse coincidence (in terms of pixels):

def pixel_coinc(binary_img1,binary_img2):

    pixel_overlap_img=binary_img1&binary_img2
    pixel_overlap_count=pixel_overlap_img.sum()
    pixel_fraction=pixel_overlap_img.sum()/binary_img1.sum()

    return pixel_overlap_img,pixel_overlap_count,pixel_fraction

# Function to analyse coincidence (in terms of features):

def feature_coinc(binary_img1,binary_img2):

    features_number,labelled_img1=label_image(binary_img1)

# Finds pixel overlap between the two images

    coinc_img=binary_img1&binary_img2

# Gives a coincident image with the pixels being equal to label

    coinc_labels=labelled_img1*coinc_img

# Counts number of unique occureences in the image

    coinc_list, coinc_pixels = np.unique(coinc_labels, return_counts=True)

    total_labels=labelled_img1.max()
    total_labels_coinc=len(coinc_list)
    fraction_coinc=total_labels_coinc/total_labels

# Fraction of overlap in each feature
# First of all, count the number of unique occurances in original image

    label_list, label_pixels = np.unique(labelled_img1, return_counts=True)
    fract_pixels_overlap=[]

    for i in range(len(coinc_list)):
        overlap_pixels=coinc_pixels[i]
        label=coinc_list[i]
        total_pixels=label_pixels[label]
        fract=1.0*overlap_pixels/total_pixels
        fract_pixels_overlap.append(fract)

# Generate the images.

    coinc_list[0]=1000000

# Generate binary image only from labels in coinc list.

    coinc_features_img=np.isin(labelled_img1,coinc_list)
    coinc_list[0]=0

# Generate image from numbers not in coinc list.

    non_coinc_features_img=~np.isin(labelled_img1,coinc_list)

    return coinc_list,coinc_pixels,fraction_coinc,coinc_features_img,non_coinc_features_img,fract_pixels_overlap

# Function to measure minimum distances between two sets of data:

def min_distance(measurements1,measurements2):

    s1 = measurements1[["centroid-0","centroid-1"]].to_numpy()
    s2 = measurements2[["centroid-0","centroid-1"]].to_numpy()
    min_lengths=spatial.distance.cdist(s1,s2).min(axis=1)

    return min_lengths

for i in range(len(pathList)):

    directory=pathList[i]+"/"

# Run functions for aptamer channel.

    filename="Apt.tif"
    apt_img=io.imread(directory+filename)
    #numpydata = asarray(apt_img)
    #apt_threshold_value=(numpydata.mean() + 3*numpydata.std())
    apt_binary_img=threshold_img_fixed(apt_img,threshold_value=19300)
    #print("Aptamer threshold value was %d."%apt_threshold_value)
    im = Image.fromarray(apt_binary_img)
    im.save(directory+'Apt_Binary.tif')
    apt_features_number,apt_labelled_img=label_image(apt_binary_img)
    print("%d features were detected in the aptamer image."%apt_features_number)
    apt_measurements=analyse_labelled_img(apt_labelled_img,apt_img)

# Run functions for antibody channel.

    filename="Ab.tif"
    ab_img=io.imread(directory+filename)
    #numpydata = asarray(ab_img)
    #ab_threshold_value=(numpydata.mean() + 3*numpydata.std())
    ab_binary_img=threshold_img_fixed(ab_img,threshold_value=13500)
    #print("Antibody threshold value was %d."%ab_threshold_value)
    im = Image.fromarray(ab_binary_img)
    im.save(directory+'Ab_Binary.tif')
    ab_features_number,ab_labelled_img=label_image(ab_binary_img)
    print("%d features were detected in the antibody image."%ab_features_number)
    ab_measurements=analyse_labelled_img(ab_labelled_img,ab_img)

# Run functions for dapi channel.

    filename="Dapi.tif"
    dapi_img=io.imread(directory+filename)
    dapi_binary_img=threshold_img_fixed(dapi_img,threshold_value=1550)
    im = Image.fromarray(dapi_binary_img)
    im.save(directory+'Dapi_Binary.tif')
    dapi_features_number,dapi_labelled_img=label_image(dapi_binary_img)
    print("%d features were detected in the dapi image."%dapi_features_number)
    dapi_measurements=analyse_labelled_img(dapi_labelled_img,dapi_img)

# Run coincidence functions.

    ab_pixel_overlap_img,ab_pixel_overlap_count,ab_pixel_fraction=pixel_coinc(ab_binary_img,apt_binary_img)
    print("%.2f of antibody pixels had coincidence with aptamer image."%ab_pixel_fraction)

    abdapi_pixel_overlap_img,abdapi_pixel_overlap_count,abdapi_pixel_fraction=pixel_coinc(ab_binary_img,dapi_binary_img)
    print("%.2f of antibody pixels had coincidence with dapi image."%abdapi_pixel_fraction)

    apt_pixel_overlap_img,apt_pixel_overlap_count,apt_pixel_fraction=pixel_coinc(apt_binary_img,ab_binary_img)
    print("%.2f of aptamer pixels had coincidence with antibody image."%apt_pixel_fraction)

    aptdapi_pixel_overlap_img,aptdapi_pixel_overlap_count,aptdapi_pixel_fraction=pixel_coinc(apt_binary_img,dapi_binary_img)
    print("%.2f of aptamer pixels had coincidence with dapi image."%aptdapi_pixel_fraction)

    apt_coinc_list,apt_coinc_pixels,apt_fraction_coinc,apt_coinc_features_img,apt_noncoinc_features_img,apt_fraction_pixels_overlap=feature_coinc(apt_binary_img,ab_binary_img)
    print("%.2f of aptamer features had coincidence with features in antibody image. Average overlap was %2f."%(apt_fraction_coinc,sum(apt_fraction_pixels_overlap)/len(apt_fraction_pixels_overlap)))

    aptdapi_coinc_list,aptdapi_coinc_pixels,aptdapi_fraction_coinc,aptdapi_coinc_features_img,aptdapi_noncoinc_features_img,aptdapi_fraction_pixels_overlap=feature_coinc(apt_binary_img,dapi_binary_img)
    print("%.2f of aptamer features had coincidence with features in dapi image. Average overlap was %2f."%(aptdapi_fraction_coinc,sum(aptdapi_fraction_pixels_overlap)/len(aptdapi_fraction_pixels_overlap)))

    abdapi_coinc_list,abdapi_coinc_pixels,abdapi_fraction_coinc,abdapi_coinc_features_img,abdapi_noncoinc_features_img,abdapi_fraction_pixels_overlap=feature_coinc(ab_binary_img,dapi_binary_img)
    print("%.2f of antibody features had coincidence with features in dapi image. Average overlap was %2f."%(abdapi_fraction_coinc,sum(abdapi_fraction_pixels_overlap)/len(abdapi_fraction_pixels_overlap)))

    apt_coinc_tosave=apt_coinc_features_img*apt_img
    im = Image.fromarray(apt_coinc_tosave)
    im.save(directory+'Aptamer_features_coinc.tif')

    apt_noncoinc_tosave=apt_noncoinc_features_img*apt_img
    im = Image.fromarray(apt_noncoinc_tosave)
    im.save(directory+'Aptamer_features_noncoinc.tif')

    ab_coinc_list,ab_coinc_pixels,ab_fraction_coinc,ab_coinc_features_img,ab_noncoinc_features_img,ab_fraction_pixels_overlap=feature_coinc(ab_binary_img,apt_binary_img)
    print("%.2f of antibody features had coincidence with features in aptamer image. Average overlap was %2f."%(ab_fraction_coinc,sum(ab_fraction_pixels_overlap)/len(ab_fraction_pixels_overlap)))

    ab_coinc_tosave=ab_coinc_features_img*ab_img
    im = Image.fromarray(ab_coinc_tosave)
    im.save(directory+'Ab_features_coinc.tif')

    ab_noncoinc_tosave=ab_noncoinc_features_img*ab_img
    im = Image.fromarray(ab_noncoinc_tosave)
    im.save(directory+'Antibody_features_noncoinc.tif')

    apt_distance_to_nuc=min_distance(apt_measurements,dapi_measurements)*pixel_size/1000.0
    ab_distance_to_nuc=min_distance(ab_measurements,dapi_measurements)*pixel_size/1000.0
    dapi_distance_to_nuc=min_distance(dapi_measurements,dapi_measurements)*pixel_size/1000.0

# Output data.

    ab_measurements['Distance to nucleus']=min_distance(ab_measurements,dapi_measurements)
    ab_measurements.to_csv(root_directory + '/' 'ab' +str(i+1)+'.csv', sep = '\t')

    apt_measurements['Distance to nucleus']=min_distance(apt_measurements,dapi_measurements)
    apt_measurements.to_csv(root_directory + '/' 'apt' +str(i+1)+'.csv', sep = '\t')

    imRGB = np.zeros((dapi_binary_img.shape[0],dapi_binary_img.shape[1],3))
    imRGB[:,:,0] = apt_binary_img
    imRGB[:,:,1] = ab_binary_img
    imRGB[:,:,2] = dapi_binary_img

    fig, ax = plt.subplots(1,1,figsize=(40, 40))
    ax.imshow(imRGB)
    plt.savefig(directory+"binary.png")

# Open file for writing:

    Out_nc=open(directory+'/'+'Aptamer_overall_stats.txt','w')
    Out_nc.write("Number of detected features = %.2f \n" %apt_features_number)
    #Out_nc.write("Threshold value = %.2f \n" %apt_threshold_value)
    Out_nc.write("Pixel coincidence = %.2f \n" %apt_pixel_fraction)
    Out_nc.write("Feature coincidence with antibody = %.2f \n" %apt_fraction_coinc)
    Out_nc.write("Feature coincidence with dapi = %.2f \n" %aptdapi_fraction_coinc)
    Out_nc.write("Feature overlap with antibody = %.2f \n" %(sum(apt_fraction_pixels_overlap)/len(apt_fraction_pixels_overlap)))
    Out_nc.write("Feature overlap with dapi = %.2f \n" %(sum(aptdapi_fraction_pixels_overlap)/len(aptdapi_fraction_pixels_overlap)))
    Out_nc.close()

    Out_nc=open(directory+'/'+'Antibody_overall_stats.txt','w')
    Out_nc.write("Number of detected features = %.2f \n" %ab_features_number)
    #Out_nc.write("Threshold value = %.2f \n" %ab_threshold_value)
    Out_nc.write("Pixel coincidence = %.2f \n" %ab_pixel_fraction)
    Out_nc.write("Feature coincidence with aptamer = %.2f \n" %ab_fraction_coinc)
    Out_nc.write("Feature coincidence with dapi = %.2f \n" %abdapi_fraction_coinc)
    Out_nc.write("Feature overlap with aptamer= %.2f \n" %(sum(ab_fraction_pixels_overlap)/len(ab_fraction_pixels_overlap)))
    Out_nc.write("Feature overlap with dapi= %.2f \n" %(sum(abdapi_fraction_pixels_overlap)/len(abdapi_fraction_pixels_overlap)))
    Out_nc.close()


    Output_overall = pd.DataFrame(columns=['Image','Number of nuclei','ab number of features','apt number of features',
                                           'ab pixel coincidence','ab feature coincidence','ab feature overlap','apt pixel coincidence','apt feature coincidence',
                                           'apt feature overlap'])

    Output_overall=Output_overall.append({'Image': directory,'Number of nuclei':dapi_features_number,'ab number of features':ab_features_number,
                                            'apt number of features':apt_features_number,'ab pixel coincidence':ab_pixel_fraction,'ab feature coincidence':ab_fraction_coinc,
                                            'ab feature overlap':(sum(ab_fraction_pixels_overlap)/len(ab_fraction_pixels_overlap)),
                                            'apt pixel coincidence':apt_pixel_fraction,'apt feature coincidence':apt_fraction_coinc,
                                            'apt feature overlap':(sum(apt_fraction_pixels_overlap)/len(apt_fraction_pixels_overlap)),
                                          }, ignore_index=True)

    Output_overall.to_csv(root_directory + '/' 'all' +str(i+1)+'.csv', sep = '\t')
