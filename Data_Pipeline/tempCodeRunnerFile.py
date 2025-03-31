image_files = []
for i in range(7):
    csv_path = os.path.join(script_dir, 'Image_Data', 'Stacked_Force_Data', f'Image_{i}_Stacked_Force.csv')
    if not os.path.isfile(csv_path):
        print(f"WARNING: File not found: {csv_path}")
    else:
        image_files.append(csv_path)