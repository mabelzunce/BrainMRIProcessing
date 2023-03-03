import csv


def extract_labels(file, start, end):
    """Read file, extract labels from START to END """

    labels = {}
    with open(file, 'r') as f:
        lines = f.readlines()

        for line in lines:
            line = line.strip().split()

            if line:  # convert to int the number of label
                try:
                    n_label = int(line[0])
                    name_label = line[1]

                    # segmentation
                    if start <= n_label <= end:
                        labels[n_label] = name_label

                except ValueError:
                    pass

    return labels




if __name__ == '__main__':
    path_fscolorlut = "/home/sol/COVID/atlas/LUT_txt/FreeSurferColorLUT.txt"
    path_asegstatslut = "/home/sol/COVID/atlas/LUT_txt/ASegStatsLUT.txt"

    segmentation = extract_labels(path_asegstatslut, 2, 255)
    parcellation_DK_atlas = extract_labels(path_fscolorlut, 1000, 2035)
    parcellation_Des_atlas = extract_labels(path_fscolorlut, 11100, 12175)

    # write csv
    path_labels_csv = "/home/sol/COVID/atlas/labels.csv"

    with open(path_labels_csv, 'w') as file:
        csv_writer = csv.writer(file)

        # headers
        headers = ['n_label', 'structure', 'atlas']
        csv_writer.writerow(headers)

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
