import os, glob
import pandas as pd

# Begin with aptamer metrics:    
# Where to save output:

root_directory=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all_apt_metrics/"

# Files to merge:

path=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all_apt_metrics/"

all_files = glob.glob(os.path.join(path, "apt*.csv"))

# Merge dfs

all_df = []
for f in all_files:
    df = pd.read_csv(f, sep=',')
    df['file'] = f.split('/')[-1]
    all_df.append(df)

merged_df = pd.concat(all_df, ignore_index=True, sort=True)

# Save merged df

merged_df.to_csv(root_directory + '/' + 'apt_merged.csv', sep = '\t')

# Now antibody metrics:
# Where to save output

root_directory=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all_ab_metrics/"

# Files to merge:

path=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all_ab_metrics/"

all_files = glob.glob(os.path.join(path, "ab*.csv"))

# Merge dfs

all_df = []
for f in all_files:
    df = pd.read_csv(f, sep=',')
    df['file'] = f.split('/')[-1]
    all_df.append(df)

merged_df = pd.concat(all_df, ignore_index=True, sort=True)

# Save merged df

merged_df.to_csv(root_directory + '/' + 'ab_merged.csv', sep = '\t')

# Now overall metrics:
# Where to save output

root_directory=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all/"

# Files to merge:

path=r"D:/ECAS_Cohort/FL_Crops/SD053_16_crops/all/"

all_files = glob.glob(os.path.join(path, "all*.csv"))

# Merge dfs

all_df = []
for f in all_files:
    df = pd.read_csv(f, sep=',')
    df['file'] = f.split('/')[-1]
    all_df.append(df)

merged_df = pd.concat(all_df, ignore_index=True, sort=True)

# Save merged df

merged_df.to_csv(root_directory + '/' + 'all_merged.csv', sep = '\t')
