import pandas as pd
import os
from sklearn.metrics import precision_score, recall_score, f1_score
import gzip
# Function to extract contig names and lengths from a .gfa file
def extract_contigs_and_lengths(gfa_path):
    contigs = []
    lengths = []
    
    # Open the file, supporting both regular and gzip-compressed files
    open_func = gzip.open if gfa_path.endswith('.gz') else open
    with open_func(gfa_path, 'rt') as gfa_file:  # 'rt' mode for reading text
        for line in gfa_file:
            if line.startswith('S'):  # 'S' lines define contigs in GFA format
                parts = line.strip().split('\t')
                contig_name = parts[1]
                contig_length = int(parts[3].split(':')[-1])  # Extract length from "LN:<length>"
                contigs.append(contig_name)
                lengths.append(contig_length)
    
    return pd.DataFrame({"Contig": contigs, "Length": lengths})
# Function to read ground truth data and keep only the 'contig' column
def read_ground_truth(gt_path):
    """
    Reads the ground truth file and creates a DataFrame where the 'GT' column 
    is 1 for plasmids and 0 for chromosomes.
    """
    # Read the file with the updated format
    ground_truth_df = pd.read_csv(gt_path, sep=',')  # CSV with comma delimiter

    # Check if required columns are present
    required_columns = {"contig", "label"}
    if not required_columns.issubset(ground_truth_df.columns):
        raise ValueError(f"Ground truth file {gt_path} does not contain the required columns: {required_columns}")

    # Normalize column names and create GT column
    ground_truth_df["Contig"] = ground_truth_df["contig"].astype(str).str.lower()  # Normalize contig names
    ground_truth_df["GT"] = ground_truth_df["label"].apply(lambda x: 1 if x == "plasmid" else 0)

    return ground_truth_df[["Contig", "GT"]]  # Return standardized DataFrame

# Function to extract numbers (X) from the Contigs column in the last file
def extract_contig_numbers(txt_path):
    contig_numbers = set()  # Use a set to handle duplicates automatically
    with open(txt_path, 'r') as file:
        for line in file:
            if not line.startswith("#"):  # Skip comment lines
                parts = line.strip().split('\t')
                if len(parts) >= 4:  # Ensure there is a Contigs column
                    contigs = parts[3]
                    for contig in contigs.split(','):
                        contig_number = contig.split(':')[0]  # Extract X from X:Y
                        contig_numbers.add(contig_number)
    return pd.DataFrame({"Contig": sorted(contig_numbers)})  # Return as a DataFrame

# Function to create the combined DataFrame for each sample
def create_combined_dataframe(contigs_df, ground_truth_df, plasbinflow_df):
    """
    Creates a combined DataFrame for the given sample, with binary indicators 
    for PlasBinFlow predictions and ground truth labels.
    """
    # Normalize contig names
    contigs_df["Contig"] = contigs_df["Contig"].astype(str).str.lower()
    ground_truth_df["Contig"] = ground_truth_df["Contig"].astype(str).str.lower()
    plasbinflow_df["Contig"] = plasbinflow_df["Contig"].astype(str).str.lower()

    # Merge contigs_df with PlasBinFlow predictions and ground truth
    contigs_df["PlasBinFlow_Pred"] = contigs_df["Contig"].apply(
        lambda x: 1 if x in plasbinflow_df["Contig"].values else 0
    )
    contigs_df["GT"] = contigs_df["Contig"].apply(
        lambda x: 1 if x in ground_truth_df[ground_truth_df["GT"] == 1]["Contig"].values else 0
    )

    return contigs_df
# Function to compute unweighted metrics
def compute_unweighted_metrics(pred_column, gt_column):
    precision = precision_score(gt_column, pred_column, zero_division=0)
    recall = recall_score(gt_column, pred_column, zero_division=0)
    f1 = f1_score(gt_column, pred_column, zero_division=0)
    return precision, recall, f1

# Function to compute weighted metrics
def compute_weighted_metrics(combined_df):
    tp_contigs = combined_df[(combined_df["PlasBinFlow_Pred"] == 1) & (combined_df["GT"] == 1)]
    fp_contigs = combined_df[(combined_df["PlasBinFlow_Pred"] == 1) & (combined_df["GT"] == 0)]
    fn_contigs = combined_df[(combined_df["PlasBinFlow_Pred"] == 0) & (combined_df["GT"] == 1)]

    tp_length = tp_contigs["Length"].sum()
    fp_length = fp_contigs["Length"].sum()
    fn_length = fn_contigs["Length"].sum()

    weighted_precision = tp_length / (tp_length + fp_length) if (tp_length + fp_length) > 0 else 0
    weighted_recall = tp_length / (tp_length + fn_length) if (tp_length + fn_length) > 0 else 0

    if weighted_precision + weighted_recall > 0:
        weighted_f1 = 2 * (weighted_precision * weighted_recall) / (weighted_precision + weighted_recall)
    else:
        weighted_f1 = 0

    return weighted_precision, weighted_recall, weighted_f1

# Function to compute granular statistics for each sample
def compute_granular_statistics(combined_df, sample_name, seed_file_path):
    total_contigs = len(combined_df)
    plasmid_contigs = combined_df[combined_df["GT"] == 1]
    predicted_contigs = combined_df[combined_df["PlasBinFlow_Pred"] == 1]

    # TP, FP, FN
    tp_contigs = predicted_contigs[predicted_contigs["GT"] == 1]
    fp_contigs = predicted_contigs[predicted_contigs["GT"] == 0]
    fn_contigs = plasmid_contigs[plasmid_contigs["PlasBinFlow_Pred"] == 0]

    # Count seeds
    if os.path.exists(seed_file_path):
        with open(seed_file_path, 'r') as seed_file:
            num_seeds = sum(1 for _ in seed_file)  # Count the lines in the seed file
    else:
        num_seeds = 0

    # Compute metrics
    unweighted_precision, unweighted_recall, unweighted_f1 = compute_unweighted_metrics(
        combined_df["PlasBinFlow_Pred"], combined_df["GT"]
    )
    weighted_precision, weighted_recall, weighted_f1 = compute_weighted_metrics(combined_df)

    return {
        "sample": sample_name,
        "#contigs": total_contigs,
        "#plasmid contigs": len(plasmid_contigs),
        "total length of plasmid contigs": plasmid_contigs["Length"].sum(),
        "#predicted plasmid contigs": len(predicted_contigs),
        "total length of predicted plasmid contigs": predicted_contigs["Length"].sum(),
        "#seed contigs": num_seeds,
        "#TP predicted plasmid contigs": len(tp_contigs),
        "total length of TP predicted plasmids": tp_contigs["Length"].sum(),
        "#FP predicted plasmid contigs": len(fp_contigs),
        "total length of FP predicted plasmid contigs": fp_contigs["Length"].sum(),
        "#FN (unpredicted) plasmid contigs": len(fn_contigs),
        "total length of FN plasmid contigs": fn_contigs["Length"].sum(),
        "unweighted precision": unweighted_precision,
        "unweighted recall": unweighted_recall,
        "unweighted F1": unweighted_f1,
        "weighted precision": weighted_precision,
        "weighted recall": weighted_recall,
        "weighted F1": weighted_f1,
    }

# Function to save granular statistics to a CSV file
def save_granular_statistics(stats_list, output_file="granular_statistics.csv"):
    df = pd.DataFrame(stats_list)
    df.to_csv(output_file, index=False)

# Function to save metrics to separate directories
def save_metrics_to_separate_dirs(sample_metrics, combined_metrics, weighted_dir="weighted_metrics_results", unweighted_dir="unweighted_metrics_results"):
    os.makedirs(weighted_dir, exist_ok=True)
    os.makedirs(unweighted_dir, exist_ok=True)

    for sample, metrics in sample_metrics.items():
        # Weighted metrics
        weighted_file_path = os.path.join(weighted_dir, f"{sample}_metrics.txt")
        with open(weighted_file_path, 'w') as file:
            file.write(f"Weighted Precision: {metrics['weighted_precision']}\n")
            file.write(f"Weighted Recall: {metrics['weighted_recall']}\n")
            file.write(f"Weighted F1: {metrics['weighted_f1']}\n")

        # Unweighted metrics
        unweighted_file_path = os.path.join(unweighted_dir, f"{sample}_metrics.txt")
        with open(unweighted_file_path, 'w') as file:
            file.write(f"Unweighted Precision: {metrics['unweighted_precision']}\n")
            file.write(f"Unweighted Recall: {metrics['unweighted_recall']}\n")
            file.write(f"Unweighted F1: {metrics['unweighted_f1']}\n")

    # Save combined metrics
    combined_weighted_file_path = os.path.join(weighted_dir, "combined_metrics.txt")
    combined_unweighted_file_path = os.path.join(unweighted_dir, "combined_metrics.txt")

    with open(combined_weighted_file_path, 'w') as file:
        file.write(f"Weighted Precision: {combined_metrics['weighted_precision']}\n")
        file.write(f"Weighted Recall: {combined_metrics['weighted_recall']}\n")
        file.write(f"Weighted F1: {combined_metrics['weighted_f1']}\n")

    with open(combined_unweighted_file_path, 'w') as file:
        file.write(f"Unweighted Precision: {combined_metrics['unweighted_precision']}\n")
        file.write(f"Unweighted Recall: {combined_metrics['unweighted_recall']}\n")
        file.write(f"Unweighted F1: {combined_metrics['unweighted_f1']}\n")

# Main function to process samples and compute granular statistics and metrics
def process_samples_with_metrics_and_statistics(file_path):
    samples_df = pd.read_csv(file_path)
    
    granular_stats = []
    combined_df_list = []  # List to aggregate all combined dataframes for "all_samples"
    sample_metrics = {}
    combined_tp_length = 0
    combined_fp_length = 0
    combined_fn_length = 0
    combined_preds = []
    combined_gts = []

    for _, row in samples_df.iterrows():
        sample = row["sample"]
        gfa_path = row["gfa_path"]
        gt_path = row["ground_truth"]
        txt_path = row["PlasBinFlow_Pred"]
        seed_file_path = row["seed_file_path"]

        if not os.path.exists(gfa_path) or not os.path.exists(gt_path) or not os.path.exists(txt_path):
            continue

        contigs_df = extract_contigs_and_lengths(gfa_path)
        ground_truth_df = read_ground_truth(gt_path)
        plasbinflow_df = extract_contig_numbers(txt_path)
        
        combined_df = create_combined_dataframe(contigs_df, ground_truth_df, plasbinflow_df)
        combined_df_list.append(combined_df)

        # Compute granular statistics
        stats = compute_granular_statistics(combined_df, sample, seed_file_path)
        granular_stats.append(stats)

        # Compute unweighted metrics
        unweighted_precision, unweighted_recall, unweighted_f1 = compute_unweighted_metrics(
            combined_df["PlasBinFlow_Pred"], combined_df["GT"]
        )

        # Compute weighted metrics
        weighted_precision, weighted_recall, weighted_f1 = compute_weighted_metrics(combined_df)

        sample_metrics[sample] = {
            "unweighted_precision": unweighted_precision,
            "unweighted_recall": unweighted_recall,
            "unweighted_f1": unweighted_f1,
            "weighted_precision": weighted_precision,
            "weighted_recall": weighted_recall,
            "weighted_f1": weighted_f1,
        }

        # Aggregate for combined metrics
        combined_preds.extend(combined_df["PlasBinFlow_Pred"].values)
        combined_gts.extend(combined_df["GT"].values)
        combined_tp_length += combined_df[(combined_df["PlasBinFlow_Pred"] == 1) & (combined_df["GT"] == 1)]["Length"].sum()
        combined_fp_length += combined_df[(combined_df["PlasBinFlow_Pred"] == 1) & (combined_df["GT"] == 0)]["Length"].sum()
        combined_fn_length += combined_df[(combined_df["PlasBinFlow_Pred"] == 0) & (combined_df["GT"] == 1)]["Length"].sum()

    # Compute combined unweighted metrics
    combined_unweighted_precision, combined_unweighted_recall, combined_unweighted_f1 = compute_unweighted_metrics(
        combined_preds, combined_gts
    )

    # Compute combined weighted metrics
    combined_weighted_precision = combined_tp_length / (combined_tp_length + combined_fp_length) if (combined_tp_length + combined_fp_length) > 0 else 0
    combined_weighted_recall = combined_tp_length / (combined_tp_length + combined_fn_length) if (combined_tp_length + combined_fn_length) > 0 else 0
    combined_weighted_f1 = 2 * (combined_weighted_precision * combined_weighted_recall) / (combined_weighted_precision + combined_weighted_recall) if (combined_weighted_precision + combined_weighted_recall) > 0 else 0

    combined_metrics = {
        "unweighted_precision": combined_unweighted_precision,
        "unweighted_recall": combined_unweighted_recall,
        "unweighted_f1": combined_unweighted_f1,
        "weighted_precision": combined_weighted_precision,
        "weighted_recall": combined_weighted_recall,
        "weighted_f1": combined_weighted_f1,
    }

    # Combine all samples' dataframes for "all_samples"
    combined_all_samples_df = pd.concat(combined_df_list, ignore_index=True)

    # Compute granular statistics for "all_samples"
    all_samples_stats = compute_granular_statistics(combined_all_samples_df, "all_samples", seed_file_path="")
    granular_stats.append(all_samples_stats)

    return granular_stats, sample_metrics, combined_metrics

# Specify the input file path
input_file_path = "/home/paa40/projects/ctb-chauvec/PLASMIDS/tools/PlasBinFlow_platon/PlasBin-flow/code/Random_50_Eval_Path.csv"

# Process the samples and compute metrics and granular statistics
granular_statistics, sample_metrics, combined_metrics = process_samples_with_metrics_and_statistics(input_file_path)

# Save metrics to separate directories
save_metrics_to_separate_dirs(sample_metrics, combined_metrics)

# Save granular statistics to a CSV file
save_granular_statistics(granular_statistics)

print("Metrics saved to 'weighted_metrics_results' and 'unweighted_metrics_results'")
print("Granular statistics saved to 'granular_statistics.csv'")

