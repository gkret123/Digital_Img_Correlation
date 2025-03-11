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