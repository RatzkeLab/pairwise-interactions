import random
import pandas as pd
from collections import defaultdict, Counter
from typing import List, Tuple, Set

# === Parameters ===
total_strains = 100
total_plates = 5
combinations_per_plate = 330
replicate_rows = 1
destination_rows = 16
destination_columns = list(range(2, 24)) # 22 columns (2–23)
combination_sizes = [2,4,6,8,12,16, 24, 48]

# === User input for source wells ===
user_input = input(f"Enter {total_strains} space-separated source wells (e.g., A1 B1 C1 ...):\n")
source_wells = user_input.strip().split()
if len(source_wells) != total_strains:
    raise ValueError(f"Expected {total_strains} wells but got {len(source_wells)}.")

# === Strain list and mapping ===
strains = [str(i + 1) for i in range(total_strains)]
strain_to_source = {strain: source_wells[i] for i, strain in enumerate(strains)}

# === Derived constants ===
destination_wells = [
    f"{chr(65 + row)}{col}" for row in range(destination_rows) for col in destination_columns
]
replicate_wells_per_plate = len(destination_columns)

# === Output containers ===
all_combinations_global: Set[Tuple[str, ...]] = set()
plate_combinations_dict = defaultdict(list)
replicate_combinations: List[Tuple[str, ...]] = []

# === Function to generate globally unique combinations with equal size distribution ===
def generate_plate_combinations(
    strains: List[str], sizes: List[int], total_needed: int, global_used: Set[Tuple[str, ...]]
) -> List[Tuple[str]]:
    combinations = []
    combinations_per_size = total_needed // len(sizes)

    for size in sizes:
        count = 0
        attempts = 0
        while count < combinations_per_size and attempts < 5000:
            comb = tuple(sorted(random.sample(strains, size)))
            if comb in global_used:
                attempts += 1
                continue
            combinations.append(comb)
            global_used.add(comb)
            count += 1
            attempts = 0 # reset on success
    return combinations

# === Generate 22 fixed replicate combinations ===
for size in combination_sizes[:replicate_wells_per_plate]:
    replicate_combinations.append(tuple(sorted(random.sample(strains, size))))

# === Generate combinations for all plates ===
for plate_id in range(total_plates):
    plate_name = f"Plate{plate_id + 1}"
    combos = generate_plate_combinations(strains, combination_sizes, combinations_per_plate, all_combinations_global)
    plate_combinations_dict[plate_name] = combos

# === Track strain usage ===
strain_usage = Counter({strain: 0 for strain in strains})
strain_volume = Counter({strain: 0.0 for strain in strains})

# === Output generation ===
for plate_name, combinations in plate_combinations_dict.items():
    picklist_data = []
    presence_absence_data = []
    combination_table_data = []
    dest_well_idx = 0

    # Unique combinations
    for row in range(destination_rows - replicate_rows):
        for col in destination_columns:
            if dest_well_idx >= len(combinations):
                break
            dest_well = destination_wells[dest_well_idx]
            comb = combinations[dest_well_idx]
            comb_str = " ".join(comb)
            volume = round(1200 / len(comb), 1)

            for strain in comb:
                picklist_data.append([strain_to_source[strain], volume, dest_well])
                strain_usage[strain] += 1
                strain_volume[strain] += volume

            combination_table_data.append([plate_name, dest_well, comb_str])
            row_data = {"Plate": plate_name, "Destination_Well": dest_well}
            for strain in strains:
                row_data[strain] = 1 if strain in comb else 0
            presence_absence_data.append(row_data)
            dest_well_idx += 1

    # Replicate combinations
    for col_idx, col in enumerate(destination_columns[:22]):
        dest_well = destination_wells[(destination_rows - replicate_rows) * len(destination_columns) + col_idx]
        comb = replicate_combinations[col_idx % len(replicate_combinations)]
        comb_str = " ".join(comb)
        volume = round(1200 / len(comb), 1)

        for strain in comb:
            picklist_data.append([strain_to_source[strain], volume, dest_well])
            strain_usage[strain] += 1
            strain_volume[strain] += volume

        combination_table_data.append([plate_name, dest_well, comb_str])
        row_data = {"Plate": plate_name, "Destination_Well": dest_well}
        for strain in strains:
            row_data[strain] = 1 if strain in comb else 0
        presence_absence_data.append(row_data)

    # Save plate-specific files
    pd.DataFrame(picklist_data, columns=["Source Well", "Transfer Volume (nL)", "Destination Well"]).to_csv(
        f"picklist_{plate_name}.csv", index=False
    )
    pd.DataFrame(presence_absence_data).to_csv(f"presence_absence_{plate_name}.csv", index=False)
    pd.DataFrame(combination_table_data, columns=["Plate", "Destination Well", "Combination"]).to_csv(
        f"combination_table_{plate_name}.csv", index=False
    )

# === Save usage summary ===
pd.DataFrame({
    "Strain_ID": list(strain_usage.keys()),
    "Total_Uses": list(strain_usage.values()),
    "Total_Volume_nL": [strain_volume[s] for s in strain_usage.keys()]
}).to_csv("strain_usage_summary.csv", index=False)

print(" All files generated successfully with equal distribution across combination sizes.")