import pandas as pd
import numpy as np
import os

# Define the input file name and the chunk size
input_file = 'Extracted Data from VIC-EDU_Poly.csv'

chunk_size = 204 # 204 lines per file becuase vic software takes a 200 data point line and adds 4 lines of header

metrics_per_line = 9 #we export 9 metrics per line from VIC-EDU software (we can add more or less if needed)
#Metrics chosen: ["X [mm]","Y [mm]","Z [mm]","U [mm]","V [mm]","W [mm]","exx [1] - Lagrange","eyy [1] - Lagrange","exy [1] - Lagrange"]

# Define the forces for each image
forces = [0, 38, 92, 140, 183, 226, np.nan] #nan is for the last image, which is after the test hasd been completed, force is zero but dont want to confuse with the first image
# Read all lines from the file
with open(input_file, 'r', encoding='utf-8') as f:
    lines = f.readlines()

total_lines = len(lines)
print(f"Total lines: {total_lines}")

# Calculate the number of chunks required, this seperated the data collected from each image
num_chunks = (total_lines + chunk_size - 1) // chunk_size

# Loop over each chunk and write it to a separate file
for i in range(num_chunks):
    # Get the chunk of 204 lines (or fewer for the last file)
    chunk_lines = lines[i * chunk_size:(i + 1) * chunk_size]
    #remove first line of each file
    chunk_lines.pop(0)
    output_file = os.path.join('Image_Data', f'Image_{i}_Data_{forces[i]}N.csv')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', encoding='utf-8', newline='') as out_f:
        out_f.writelines(chunk_lines)
    print(f"Wrote {output_file} with {len(chunk_lines)} lines.")

#the first line of each file should be the keys, the next line is the header, however, note that the header repeats itself for each key. work with the inidividual files to prep the dataset
#read the first file to get the keys
df_extract_headers = pd.read_csv(os.path.join('Image_Data', 'Image_0_Data_0N.csv'), index_col=0)
first_level_header = [f"Image{index}" for index in range(num_chunks)]
#extract the keys and drop unnamed values
second_level_header = [k for k in df_extract_headers.keys() if 'Unnamed' not in k]
print(second_level_header)
#extract the third level headers
third_level_header = df_extract_headers.iloc[0].unique()
print(third_level_header)
df = pd.DataFrame(columns = pd.MultiIndex.from_tuples([(H1, H2, H3) for H1, H2, H3 in zip(first_level_header, second_level_header, third_level_header)]))




#make mini dfs for each L1,L2,L3,etc. in each file, then concatonate them to the main df
df_image_0 = pd.read_csv(os.path.join('Image_Data', 'Image_0_Data_0N.csv'), index_col=0)
#df_image_0 = df_image_0.drop(df_image_0.index[0])
ff_image0_groups = {}
num_groups = df_image_0.shape[1] // metrics_per_line

ff_all_images = {}
#loop over all the files and give each file a key from second level header
for img_idx in range(num_chunks):
    df_image = pd.read_csv(f'Image_{img_idx}_Data_{forces[img_idx]}N.csv', index_col=0)
    num_groups = df_image.shape[1] // metrics_per_line
    image_groups = {}
    #loop over all the groups in the file
    for grp in range(num_groups):
        group = df_image.iloc[:, grp * metrics_per_line:(grp + 1) * metrics_per_line]
        group.columns = third_level_header
        group = group.drop(group.index[0]).reset_index(drop=True)
        image_groups[f"ff_image_L{grp}"] = group
    ff_all_images[f"Image_{img_idx}"] = image_groups

print(ff_all_images["Image_3"]["ff_image_L1"])

#add a 2 column to each group, force and L value
for img_idx in range(num_chunks):
    for grp in range(num_groups):
        ff_all_images[f"Image_{img_idx}"][f"ff_image_L{grp}"]['Force'] = forces[img_idx]
        ff_all_images[f"Image_{img_idx}"][f"ff_image_L{grp}"]['Line'] = grp
        
print(ff_all_images["Image_1"]["ff_image_L1"])

#stack all the groups of the same force value into a single df
ff_all_images_stacked_Force = {}
for img_idx in range(num_chunks):
    stacked_df = pd.concat([ff_all_images[f"Image_{img_idx}"][f"ff_image_L{grp}"] for grp in range(num_groups)])
    ff_all_images_stacked_Force[f"Image_{img_idx}"] = stacked_df

#stack all the groups with the same L value into a single df
ff_all_images_stacked_L = {}
for grp in range(num_groups):
    stacked_df_L = pd.concat([ff_all_images[f"Image_{img_idx}"][f"ff_image_L{grp}"] for img_idx in range(num_chunks)])
    ff_all_images_stacked_L[f"ff_image_L{grp}"] = stacked_df_L
#print(ff_all_images_stacked_L["ff_image_L1"])
