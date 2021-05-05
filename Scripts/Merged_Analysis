import pandas as pd
import matplotlib.pyplot as plt

directory=r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all_ab_metrics/"

# Read merged antibody df from csv file

df = pd.read_csv(r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all_ab_metrics/ab_merged.csv")

data = df['Distance to nucleus'], df['area'], df['major_axis_length'], 
df['minor_axis_length'], df['perimeter'], df['mean_intensity'], df['max_intensity']

# Get antibody stats

mean_ab_nuc=df['Distance to nucleus'].mean()
std_ab_nuc=df['Distance to nucleus'].std()
mean_ab_area = df['area'].mean()
std_ab_area = df['area'].std()
mean_ab_major_length = df['major_axis_length'].mean()
std_ab_major_length = df['major_axis_length'].std()
mean_ab_minor_length = df['minor_axis_length'].mean()
std_ab_minor_length = df['minor_axis_length'].std()
mean_ab_perimeter=df['perimeter'].mean()
std_ab_perimeter=df['perimeter'].std()
mean_ab_intensities=df['mean_intensity'].mean()
std_ab_intensities=df['mean_intensity'].std()
mean_ab_max_intensities=df['max_intensity'].mean()
std_ab_max_intensities=df['max_intensity'].std()

# Record antibody stats in a txt file

Out_nc=open(directory+'/'+'ab_merged_overall_stats.txt','w')
Out_nc.write("Area = %.2f +/- %.2f \n" %(mean_ab_area,std_ab_area))
Out_nc.write("Major Axis Length = %.2f +/- %.2f \n" %(mean_ab_major_length,std_ab_major_length))
Out_nc.write("Minor Axis Length = %.2f +/- %.2f \n" %(mean_ab_minor_length,std_ab_minor_length))
Out_nc.write("Perimeter = %.2f +/- %.2f \n" %(mean_ab_perimeter,std_ab_perimeter))
Out_nc.write("Mean intensity = %.2f +/- %.2f \n" %(mean_ab_intensities,std_ab_intensities))
Out_nc.write("Max intensity = %.2f +/- %.2f \n" %(mean_ab_max_intensities,std_ab_max_intensities))
Out_nc.write("Closest nucleus = %.2f +/- %.2f \n" %(mean_ab_nuc,std_ab_nuc))
Out_nc.close()

# Plot some antibody histograms

plt.hist(df['Distance to nucleus'], bins = 100,range=[0,100], rwidth=0.9,color='green')
plt.xlabel('Distance to Nearest Nucleus (\u03bcm)')
plt.ylabel('Number of Features')
plt.title('Antibody Distance to Nearest Nucleus')
#plt.ylim([0,400])
plt.savefig(directory+'/'+'ab_distance_to_nuc.pdf')
plt.show()

plt.hist(df['area'], bins = 100,range=[5,100], rwidth=0.9,color='green')
plt.xlabel('Area of feature (\u03bcm$^2$)')
plt.ylabel('Number of Features')
plt.title('Area of Antibody Feature')
plt.ylim([0,50])
plt.savefig(directory+'/'+'antibody_areas.pdf')
plt.show()

plt.hist(df['major_axis_length'], bins = 100,range=[1,100], rwidth=0.9,color='green')
plt.xlabel('Length (major axis) (\u03bcm)')
plt.ylabel('Number of Features')
plt.title('Length of Antibody Features')
#plt.ylim([0,400])
plt.savefig(directory+'/'+'antibody_lengths.pdf')
plt.show()

plt.hist(df['perimeter'], bins = 100,range=[1,100], rwidth=0.9,color='green')
plt.xlabel('Perimeter of Feature (\u03bcm)')
plt.ylabel('Number of Features')
plt.title('Perimeter of Antibody Features')
#plt.ylim([0,400])
plt.savefig(directory+'/'+'antibody_perimeters.pdf')
plt.show()

plt.hist(df['mean_intensity'], bins = 100,range=[13500,30000], rwidth=0.9,color='green')
plt.xlabel('Mean Intensity (AU)')
plt.ylabel('Number of Features')
plt.title('Antibody Mean Intensities')
#plt.ylim([0,67])
plt.savefig(directory+'/'+'antibody_intensities_mean.pdf')
plt.show()

plt.hist(df['max_intensity'], bins = 100,range=[13500,30000], rwidth=0.9,color='green')
plt.xlabel('Maximum Intensity (AU)')
plt.ylabel('Number of Features')
plt.title('Antibody Maximum Intensities')
#plt.ylim([0,75])
plt.savefig(directory+'/'+'antibody_intensities_max.pdf')
plt.show()

# Now do the same for aptamer df

directory=r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all_apt_metrics/"

df = pd.read_csv(r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all_apt_metrics/apt_merged.csv")

data=df['Distance to nucleus'], df['area'], df['major_axis_length'], df['minor_axis_length'],
df['perimeter'], df['mean_intensity'], df['max_intensity']

# Get some aptamer stats

mean_apt_nuc=df['Distance to nucleus'].mean()
std_apt_nuc=df['Distance to nucleus'].std()
mean_apt_area = df['area'].mean()
std_apt_area = df['area'].std()
mean_apt_major_length = df['major_axis_length'].mean()
std_apt_major_length = df['major_axis_length'].std()
mean_apt_minor_length = df['minor_axis_length'].mean()
std_apt_minor_length = df['minor_axis_length'].std()
mean_apt_perimeter=df['perimeter'].mean()
std_apt_perimeter=df['perimeter'].std()
mean_apt_intensities=df['mean_intensity'].mean()
std_apt_intensities=df['mean_intensity'].std()
mean_apt_max_intensities=df['max_intensity'].mean()
std_apt_max_intensities=df['max_intensity'].std()

# Record aptamer stats in a txt file

Out_nc=open(directory+'/'+'aptamer_merged_overall_stats.txt','w')
Out_nc.write("Area = %.2f +/- %.2f \n" %(mean_ab_area,std_ab_area))
Out_nc.write("Major Axis Length = %.2f +/- %.2f \n" %(mean_ab_major_length,std_ab_major_length))
Out_nc.write("Minor Axis Length = %.2f +/- %.2f \n" %(mean_ab_minor_length,std_ab_minor_length))
Out_nc.write("Perimeter = %.2f +/- %.2f \n" %(mean_ab_perimeter,std_ab_perimeter))
Out_nc.write("Mean intensity = %.2f +/- %.2f \n" %(mean_ab_intensities,std_ab_intensities))
Out_nc.write("Max intensity = %.2f +/- %.2f \n" %(mean_ab_max_intensities,std_ab_max_intensities))
Out_nc.write("Closest nucleus = %.2f +/- %.2f \n" %(mean_ab_nuc,std_ab_nuc))
Out_nc.close()

# Plot some aptamer histograms

plt.hist(df['Distance to nucleus'], bins = 100,range=[1,100], rwidth=0.9,color='red')
plt.xlabel('Distance to Nearest Nucleus (\u03bcm)')
plt.ylabel('Number of Features')
#plt.ylim([0,400])
plt.title('Aptamer Distance to Nearest Nucleus')
plt.savefig(directory+'/'+'apt_distance_to_nuc.pdf')
plt.show()

plt.hist(df['area'], bins = 100,range=[5,100], rwidth=0.9,color='red')
plt.xlabel('Area of feature (\u03bcm$^2$)')
plt.ylabel('Number of Features')
plt.ylim([0,50])
plt.title('Area of Aptamer Feature')
plt.savefig(directory+'/'+'apt_areas.pdf')
plt.show()

plt.hist(df['major_axis_length'], bins = 100,range=[1,100], rwidth=0.9,color='red')
plt.xlabel('Length (major axis) (\u03bcm)')
plt.ylabel('Number of Features')
#plt.ylim([0,400])
plt.title('Length of Aptamer Features')
plt.savefig(directory+'/'+'apt_lengths.pdf')
plt.show()

plt.hist(df['perimeter'], bins = 100,range=[1,100], rwidth=0.9,color='red')
plt.xlabel('Perimeter of Feature (\u03bcm)')
plt.ylabel('Number of Features')
#plt.ylim([0,400])
plt.title('Perimeter of Aptamer Features')
plt.savefig(directory+'/'+'apt_perimeters.pdf')
plt.show()

plt.hist(df['mean_intensity'], bins = 100,range=[19300,30000], rwidth=0.9,color='red')
plt.xlabel('Mean Intensity (AU)')
plt.ylabel('Number of Features')
#plt.ylim([0,400])
plt.title('Aptamer Mean Intensities')
plt.savefig(directory+'/'+'apt_intensities_mean.pdf')
plt.show()

plt.hist(df['max_intensity'], bins = 100,range=[19300,30000], rwidth=0.9,color='red')
plt.xlabel('Maximum Intensity (AU)')
plt.ylabel('Number of Features')
#plt.ylim([0,225])
plt.title('Aptamer Maximum Intensities')
plt.savefig(directory+'/'+'apt_intensities_max.pdf')
plt.show()

# Now get some stats for overall sample data

directory=r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all"

df = pd.read_csv(r"D:/ECAS_Cohort/FL_Crops/None/SD059_14_crops/all/all_merged.csv")
df1=df.dropna()

data=df1['Number of nuclei'], df1['ab number of features'], df1['apt number of features'],
df1['ab pixel coincidence'], df1['ab feature coincidence'],
df1['ab feature overlap'], df1['apt pixel coincidence'], df1['apt feature coincidence'], df1['apt feature overlap']

mean_number_of_nuc = df['Number of nuclei'].mean()
std_number_of_nuc = df['Number of nuclei'].std()
total_number_of_nuc = df['Number of nuclei'].sum()
mean_ab_features = df['ab number of features'].mean()
std_ab_features = df['ab number of features'].std()
total_ab_features = df['ab number of features'].sum()
mean_apt_features = df['apt number of features'].mean()
std_apt_features = df['apt number of features'].std()
total_apt_features = df['apt number of features'].sum()
mean_ab_pixel_coinc = df['ab pixel coincidence'].mean()
std_ab_pixel_coinc = df['ab pixel coincidence'].std()
mean_ab_feature_coinc = df['ab feature coincidence'].mean()
std_ab_feature_coinc = df['ab feature coincidence'].std()
mean_ab_feature_overlap = df['ab feature overlap'].mean()
std_ab_feature_overlap = df['ab feature overlap'].std()
mean_apt_pixel_coinc = df['apt pixel coincidence'].mean()
std_apt_pixel_coinc = df['apt pixel coincidence'].std()
mean_apt_feature_coinc = df['apt feature coincidence'].mean()
std_apt_feature_coinc = df['apt feature coincidence'].std()
mean_apt_feature_overlap = df['apt feature overlap'].mean()
std_apt_feature_overlap = df['apt feature overlap'].std()

# Save overall stats in a txt file

Out_nc=open(directory+'/'+'merged_overall_stats.txt','w')
Out_nc.write("Total number of nuclei = %.2f \n" %(total_number_of_nuc))
Out_nc.write("Mean number of nuclei = %.2f +/- %.2f \n" %(mean_number_of_nuc,std_number_of_nuc))
Out_nc.write("Total number of ab features = %.2f \n" %(total_ab_features))
Out_nc.write("Mean number of ab features = %.2f +/- %.2f \n" %(mean_ab_features,std_ab_features))
Out_nc.write("Total number of apt features = %.2f \n" %(total_apt_features))
Out_nc.write("Mean number of apt features = %.2f +/- %.2f \n" %(mean_apt_features,std_apt_features))
Out_nc.write("ab pixel coincidence = %.2f +/- %.2f \n" %(mean_ab_pixel_coinc,std_ab_pixel_coinc))
Out_nc.write("ab feature coincidence = %.2f +/- %.2f \n" %(mean_ab_feature_coinc,std_ab_feature_coinc))
Out_nc.write("ab feature overlap = %.2f +/- %.2f \n" %(mean_ab_feature_overlap,std_ab_feature_overlap))
Out_nc.write("apt pixel coincidence = %.2f +/- %.2f \n" %(mean_apt_pixel_coinc,std_apt_pixel_coinc))
Out_nc.write("apt feature coincidence = %.2f +/- %.2f \n" %(mean_apt_feature_coinc,std_apt_feature_coinc))
Out_nc.write("apt feature overlap = %.2f +/- %.2f \n" %(mean_apt_feature_overlap,std_apt_feature_overlap))
Out_nc.close()
