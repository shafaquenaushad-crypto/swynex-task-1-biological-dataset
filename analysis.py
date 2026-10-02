import pandas as pd
from scipy.io import mmread
import gzip
import glob

print("Starting Task 2 analysis...")

# Find the dataset files
matrix_files = glob.glob("*matrix.mtx.gz")
features_files = glob.glob("*features.tsv.gz")
barcodes_files = glob.glob("*barcodes.tsv.gz")

print("\nFiles found:")
print("Matrix:", matrix_files)
print("Features:", features_files)
print("Barcodes:", barcodes_files)

# Check files
if not matrix_files:
    raise FileNotFoundError("Matrix file was not found.")

if not features_files:
    raise FileNotFoundError("Features file was not found.")

if not barcodes_files:
    raise FileNotFoundError("Barcodes file was not found.")

# Select files
matrix_file = matrix_files[0]
features_file = features_files[0]
barcodes_file = barcodes_files[0]

print("\nSelected files:")
print(matrix_file)
print(features_file)
print(barcodes_file)

# Read matrix
print("\nReading expression matrix...")

with gzip.open(matrix_file, "rb") as f:
    matrix = mmread(f)

print("Matrix loaded successfully.")

# Read features
print("Reading features...")

features = pd.read_csv(
    features_file,
    sep="\t",
    header=None,
    compression="gzip"
)

print("Features loaded successfully.")

# Read barcodes
print("Reading barcodes...")

barcodes = pd.read_csv(
    barcodes_file,
    sep="\t",
    header=None,
    compression="gzip"
)

print("Barcodes loaded successfully.")

# Basic information
print("\n========== BASIC DATA INFORMATION ==========")

print("Matrix shape:", matrix.shape)

print("Number of features/genes:", len(features))

print("Number of barcodes/spots:", len(barcodes))

# Dimension check
print("\n========== DIMENSION CHECK ==========")

print("Matrix rows:", matrix.shape[0])
print("Matrix columns:", matrix.shape[1])

if matrix.shape[0] == len(features):
    print("Rows match the number of features.")
else:
    print("Warning: Matrix rows do not match features.")

if matrix.shape[1] == len(barcodes):
    print("Columns match the number of barcodes.")
else:
    print("Warning: Matrix columns do not match barcodes.")

# Matrix statistics
print("\n========== MATRIX STATISTICS ==========")

print("Total number of entries:", matrix.shape[0] * matrix.shape[1])

print("Number of non-zero entries:", matrix.nnz)

print("Total expression/counts:", matrix.sum())

print("Minimum value:", matrix.min())

print("Maximum value:", matrix.max())

# Counts per barcode
print("\n========== BARCODE/SPOT STATISTICS ==========")

counts_per_barcode = matrix.sum(axis=0).A1

print("Minimum counts per barcode:", counts_per_barcode.min())

print("Maximum counts per barcode:", counts_per_barcode.max())

print("Mean counts per barcode:", counts_per_barcode.mean())

print(
    "Median counts per barcode:",
    pd.Series(counts_per_barcode).median()
)

# Detected genes per barcode
genes_per_barcode = (matrix > 0).sum(axis=0).A1

print("\n========== DETECTED GENES ==========")

print("Minimum detected genes:", genes_per_barcode.min())

print("Maximum detected genes:", genes_per_barcode.max())

print("Mean detected genes:", genes_per_barcode.mean())

print(
    "Median detected genes:",
    pd.Series(genes_per_barcode).median()
)

# Quality summary
quality_summary = pd.DataFrame({
    "Barcode": barcodes.iloc[:, 0],
    "Total_Counts": counts_per_barcode,
    "Detected_Genes": genes_per_barcode
})

print("\n========== QUALITY SUMMARY ==========")

print(quality_summary.head())

# Save result
quality_summary.to_csv(
    "quality_summary.csv",
    index=False
)

print("\nQuality summary saved as: quality_summary.csv")

print("\n==========================================")
print("TASK 2 BASIC ANALYSIS COMPLETED")
print("==========================================")
print()
import matplotlib.pyplot as plt

# ==========================================
# GRAPHS FOR BASIC QUALITY ANALYSIS
# ==========================================

# Graph 1: Total counts per barcode
plt.figure(figsize=(10, 5))
plt.hist(counts_per_barcode, bins=50)
plt.xlabel("Total Counts per Barcode")
plt.ylabel("Number of Barcodes")
plt.title("Distribution of Total Counts per Barcode")
plt.tight_layout()
plt.savefig("total_counts_distribution.png")
plt.show()

# Graph 2: Detected genes per barcode
plt.figure(figsize=(10, 5))
plt.hist(genes_per_barcode, bins=50)
plt.xlabel("Detected Genes per Barcode")
plt.ylabel("Number of Barcodes")
plt.title("Distribution of Detected Genes per Barcode")
plt.tight_layout()
plt.savefig("detected_genes_distribution.png")
plt.show()

print("\nGraphs saved successfully:")
print("1. total_counts_distribution.png")
print("2. detected_genes_distribution.png")