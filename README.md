# Digital_Image_Correlation

# Overview
This project processes and analyzes raw data collected from VIC-EDU experiments. The preprocessing script converts a single, large CSV file into organized, per-image CSV files. These files are structured by image (force level) and line (measurement group) to support further analysis in later stages of the pipeline.

# Notebooks Overview

1. **Preprocessing_VIC_Data.py** - Data Preprocessing:
    - Input: Raw data from VIC-EDU experiments
    - Output: Preprocessed data in separate dataframes for each image and line
    - This script preprocesses the raw data from VIC-EDU experiments. 
    - The script reads in the raw data from a CSV file and organizes it into separate CSV files for each image (force level) and line (measurement group).

2. **Data_Analysis_VIC_Data.py** - Data Analysis (In Progress): 
    - Input: Preprocessed data from VIC-EDU experiments
    - Output: Data analysis and visualization
    - This notebook analyzes the preprocessed data from VIC-EDU experiments. 
    - The notebook reads in the preprocessed data from separate CSV files for each image and line and performs various analyses and visualizations to explore the data.

3. **Material_Analysis_VIC_Data.py** - Material Analysis (In Progress):
    - Input: Preprocessed data from VIC-EDU experiments, Results from Data Analysis, and Material Properties
    - Output: Material analysis and visualizations
    - This notebook analyzes the preprocessed data from VIC-EDU experiments to extract material properties. 
    - The notebook reads in the preprocessed data from separate CSV files for each image and line and performs various analyses to extract material properties.
    - TBD: Potential analyses include stress-strain curves, Young's modulus, Poisson's ratio, ultimate tensile strength, strain heatmaps, etc.


# Setup Instructions

To run the notebooks, ensure you have the required libraries installed.

This project relies on the following libraries and tools:

1. **Requirements and Libraries:** 
    - Numpy
    - Pandas
    - Matplotlib
    - Seaborn
    - OS


2. **Tools:**
    - Python 3.x


# Data Source

VIC-3D is a non-contact, optical measurement system that uses digital image correlation (DIC) to measure full-field displacements and strains on the surface of an object under load. The system uses a stereo pair of cameras to capture images of the object from two different angles. The images are then processed to determine the 3D coordinates of points on the object's surface. The system can be used to measure the deformation of a wide range of materials, including metals, plastics, composites, and biological tissues. See Documentation Folder for more details.

# Authors and Contact

For further questions or support, contact:   
    - Gabriel Kret [gabriel.kret@cooper.edu]  
    - Maria Alvarado [maria.alvarado@cooper.edu]  
    - Jonas Margono [jonas.margono@cooper.edu]  

# File Structure

TBD
