import csv

path_fscolorlut = "/home/sol/COVID/atlas/LUT_txt/FreeSurferColorLUT.txt"
path_asegstatslut = "/home/sol/COVID/atlas/LUT_txt/ASegStatsLUT.txt"

segmentation = {}
parcellation_DK_atlas = {}
parcellation_Des_atlas = {}

# segmentation
with open(path_asegstatslut, 'r') as f:
    lines = f.readlines()

    for line in lines:
        line = line.strip().split()

        if line:  # convert to int the number of label
            try:
                n_label = int(line[0])
                name_label = line[1]

                # segmentation
                if n_label >= 2 and n_label <= 255:
                    segmentation[n_label] = name_label

            except ValueError:
                pass


# parcellation
with open(path_fscolorlut, 'r') as f:
    lines = f.readlines()

    for line in lines:
        line = line.strip().split()

        if line: # convert to int the number of label
            try:
                n_label = int(line[0])
                name_label = line[1]

                # parcellation DK atlas
                if n_label >= 1000 and n_label <= 1035:
                    parcellation_DK_atlas[n_label] = name_label

                # parcellation Des atlas
                if n_label >= 11100 and n_label <= 12175:
                    parcellation_Des_atlas[n_label] = name_label

            except ValueError:
                pass
# write csv
path_labels_csv = "/home/sol/COVID/atlas/labels.csv"

with open(path_labels_csv, 'w') as file:
    csv_writer = csv.writer(file)

    # segmentation rows
    for label in segmentation:
        row = [label, segmentation[label], 'segmentation']
        csv_writer.writerow(row)

    # parcellation DK rows
    for label in parcellation_DK_atlas:
        row = [label, parcellation_DK_atlas[label], 'parcellation_DK_atlas']
        csv_writer.writerow(row)

    # parcellation Des rows
    for label in parcellation_Des_atlas:
        row = [label, parcellation_Des_atlas[label], 'parcellation_Des_atlas']
        csv_writer.writerow(row)