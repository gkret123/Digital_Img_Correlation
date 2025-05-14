# Digital Image Correlation (DIC) – 3-Point Bending Analysis

## Overview

This repository contains the complete experimental pipeline for analyzing the linear elastic behavior of aluminum beams under three-point bending using Digital Image Correlation (DIC), analytical methods, and finite element simulations. This work was conducted as part of ME-360 at The Cooper Union.

The project integrates:
- Full-field strain measurements from the VIC-3D DIC system
- Analytical modeling using Python
- MATLAB-based postprocessing and heatmap generation
- Finite Element Analysis (FEA) using SolidWorks

## Repository Structure
├── Data_Pipeline
|     ├── DIC_CrossSection_Analysis
|     |     ├── DIC Data.zip
|     |     ├── HollowBox_Compare_Exx_Exy.py
|     |     ├── LargeBox_Compare_Exx_Exy.py
|     |     ├── SmallBox_Compare_Exx_Exy.py
|     ├── DIC_DataProcessing_MATLAB
|     |     ├── DIC_DataProcessing.m
|     |     ├── HollowBox
|     |     |     ├── DIC_DataProcessing_HollowBoxAluminum.m
|     |     |     ├── DIC_HollowBox_1000N.xlsx
|     |     |     ├── DIC_HollowBox_4000N.xlsx
|     |     ├── LargeBox
|     |     |     ├── DIC_DataProcessing_LargeBoxAluminum.m
|     |     |     ├── DIC_LargeBox_1000N.xlsx
|     |     |     ├── DIC_LargeBox_4000N.xlsx
|     |     ├── SmallBox
|     |     |     ├── DIC_DataProcessing_SmallBoxAluminum.m
|     |     |     ├── DIC_SmallBox_1000N.xlsx
|     |     |     ├── DIC_SmallBox_4000N.xlsx
|     ├── MechMat_Soln_3ptBending.py
|     ├── Python Line Data Extraction
|     |     ├── Data_Analysis_VIC_Data.py
|     |     ├── Image_Data
|     |     |     ├── Image_0_Data_0N.csv
|     |     |     ├── Image_1_Data_38N.csv
|     |     |     ├── Image_2_Data_92N.csv
|     |     |     ├── Image_3_Data_140N.csv
|     |     |     ├── Image_4_Data_183N.csv
|     |     |     ├── Image_5_Data_226N.csv
|     |     |     ├── Image_6_Data_nanN.csv
|     |     |     ├── Stacked_Force_Data
|     |     |     |     ├── Image_0_Stacked_Force.csv
|     |     |     |     ├── Image_1_Stacked_Force.csv
|     |     |     |     ├── Image_2_Stacked_Force.csv
|     |     |     |     ├── Image_3_Stacked_Force.csv
|     |     |     |     ├── Image_4_Stacked_Force.csv
|     |     |     |     ├── Image_5_Stacked_Force.csv
|     |     |     |     ├── Image_6_Stacked_Force.csv
|     |     |     |     ├── ff_all_images_stacked_Force_keys.json
|     |     |     ├── Stacked_line_Data
|     |     |     |     ├── ff_all_images_stacked_L_keys.json
|     |     |     |     ├── ff_image_L0_Stacked_Force.csv
|     |     |     |     ├── ff_image_L1_Stacked_Force.csv
|     |     |     |     ├── ff_image_L2_Stacked_Force.csv
|     |     |     |     ├── ff_image_L3_Stacked_Force.csv
|     |     |     |     ├── ff_image_L4_Stacked_Force.csv
|     |     |     |     ├── ff_image_L5_Stacked_Force.csv
|     |     |     |     ├── ff_image_L6_Stacked_Force.csv
|     |     |     |     ├── ff_image_L7_Stacked_Force.csv
|     |     |     |     ├── ff_image_L8_Stacked_Force.csv
|     |     |     |     ├── ff_image_L9_Stacked_Force.csv
|     |     ├── Input_Data
|     |     |     ├── Extracted Data from VIC-EDU_Poly.csv
|     |     ├── Plots
|     |     |     ├── _heatmap_strain_distributions_entire_beam.png
|     |     |     ├── deflection_curve.png
|     |     |     ├── shear_bending_diagrams.png
|     |     |     ├── strain_distributions_x=0.0300 m.png
|     |     |     ├── stress_distributions_x=0.0300 m.png
|     |     ├── Preprocessing_VIC_Data.py
|     ├── __pycache__
|     |     ├── MechMat_Plot_Helper_Funcs.cpython-312.pyc
|     |     ├── MechMat_Soln_3ptBending.cpython-312.pyc
|     |     ├── Preprocessing_VIC_Data.cpython-312.pyc
|     ├── tempCodeRunnerFile.py
├── Documentation
|     ├── A2_Alvarado_Kret_Margono_BOM_DIC.xlsx
|     ├── A2_M3_DIC_ReportFinal_Alvardo_Kret_Margono.pdf
|     ├── A2_M3_FinalPresentation_DIC.pptx
|     ├── A2_M3_Report_Final_Alvardo_Kret_Margono.docx
|     ├── Flowdiagram Measurements.pptx
|     ├── Kret_Alvarado_Margono_DIC_M3_Project_Definition_Updated_2.3.2025.pptx
|     ├── Materials Research.xlsx
|     ├── VIC-EDU_2022.pdf
|     ├── VIC-EDU_manual_2024_digital.pdf
├── LICENSE.docs
├── LICENSE.md
├── README.md


## Experimental Setup

- **Material**: 6061 Aluminum
- **Beam Types**:
  - Small Solid Beam
  - Hollow Box Beam
  - Large Solid Beam
- **Test Method**: 3-point bending using Instron 34TM-50
- **Strain Capture**: Stereo DIC (VIC-EDU by Correlated Solutions)
- **Loading**: Up to 4000 N, with focus on linear elastic regime

## Analysis Pipeline

### DIC Data Processing

- Images captured at 20 load increments
- Exported strain fields from VIC-3D
- MATLAB script (`DIC_DataProcessing.m`) filters data and generates:
  - Axial strain (εxx)
  - Lateral strain (εyy)
  - Shear strain (εxy)
  - Heatmaps and contour plots

### Mechanics of Materials (Python)

- `MechMat_Soln_3ptBending.py` implements beam theory for a centrally loaded simply supported beam
- Computes:
  - Bending moment and shear force diagrams
  - Vertical deflection
  - Cross-sectional strain and stress fields
- Includes yield analysis to verify linear behavior

### Finite Element Analysis (SolidWorks)

- Static simulation of the same beam and loading configuration
- Validated against both DIC and analytical models
- Visualizations include:
  - Full-length strain fields
  - Deformed and undeformed mesh comparisons

### Python Comparison Scripts

Each beam has a dedicated Python script (e.g., `HollowBox_Compare_Exx_Exy.py`) that:
- Loads DIC strain data
- Extracts cross-sectional slices at a fixed x-location
- Fits polynomial curves to the strain profile
- Plots raw data, fits, and residuals

## Results Summary

- Strong agreement in εxx across DIC, theory, and simulation
- εyy and εxy generally small in magnitude, with higher noise
- DIC performance is best on beams with more deflection (e.g., small beam at higher loads)
- High-resolution full-field strain measurements validate beam theory within expected tolerances

## Requirements

### Python

- Python 3.9+
- Dependencies:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `scipy`

### MATLAB

- MATLAB R2021a or newer
- Requires `readtable` and `scatter` plotting functions

### Other Tools

- VIC-Snap and VIC-3D from Correlated Solutions
- SolidWorks with Simulation module
- Instron testing machine

## License

- Code: [MIT License](./LICENSE)
- Data and documentation: [CC BY 4.0](./LICENSE.content)

## Authors

- Gabriel Kret – gabriel.kret@cooper.edu  
- Maria Alvarado – maria.alvarado@cooper.edu  
- Jonas Margono – jonas.margono@cooper.edu  

## Reference

For a complete explanation of the experiment, methodology, results, and discussion, see the final report in:

`Documentation/A2_M3_DIC_ReportFinal.pdf`
