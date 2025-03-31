"""
 Data Analysis Script for Module 3 Project
Team A2: Digital Image Correlation
Authors: Gabriel Kret, Maria Alvarado, Jonas Margono

Course: ME-360-A
Instructors: Dr. David Wootton, Dr. Kamau Wright
Date: Spring 2025

Explanation:
"""

import pandas as pd
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import seaborn as sns
import Preprocessing_VIC_Data as PP

script_dir = os.path.dirname(os.path.abspath(__file__))
"""
for key, value in PP.ff_all_images.items():
    print(f"Key: {key}, Value: {value}")

#export each key, data pair to a csv file labeling the file with the key (image number) and the line number (L_)

image0_L0 = PP.ff_all_images.get("Image_0").get("ff_image_L0")
print("IMAGE 0 L0:", image0_L0)"""

#find the maximum and minimum exx and eyy for each image csv in the Stacked_Force_Data folder
image_files = []
for i in range(7):
    csv_path = os.path.join(script_dir, 'Image_Data', 'Stacked_Force_Data', f'Image_{i}_Stacked_Force.csv')
    if not os.path.isfile(csv_path):
        print(f"WARNING: File not found: {csv_path}")
    else:
        image_files.append(csv_path)

def get_max_min_exx_eyy(image_file):
    df = pd.read_csv(image_file)
    max_exx = df['exx [1] - Lagrange'].max()
    min_exx = df['exx [1] - Lagrange'].min()
    max_eyy = df['eyy [1] - Lagrange'].max()
    min_eyy = df['eyy [1] - Lagrange'].min()
    return max_exx, min_exx, max_eyy, min_eyy

for index, image_file in enumerate(image_files):
    max_exx, min_exx, max_eyy, min_eyy = get_max_min_exx_eyy(image_file)
    print(f"Image: {index}")
    print(f"Max exx: {max_exx}, Min exx: {min_exx}, Max eyy: {max_eyy}, Min eyy: {min_eyy}\n")

