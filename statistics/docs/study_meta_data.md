# Study meta data

### *!!! IMPORTANT: MAJOR updates made to CellTypes format as of 2022-10-08 !!!*

******

#### General file rules

* all file names should start with prefix [StudyName]
* all files must contain at minimum the columns described below but may contain additional columns as well
* a column must exist for any clinical variables to be analyzed, including variables to be compared and variables on which data will be filtered
* while not a requirement, it is strongly recommended that column names NOT contain any spaces or special characters
* column names may NOT start with a number
* repeated values within a column must match EXACTLY including case and spaces
* no values should contain leading or trailing whitespace
* column names should not be used in multiple XLSX files except for variables to used for joining tables (CellDive_ID, Patient_ID)

******

#### [StudyName]_Samples.xlsx
All meta data corresponding to individual samples
| Column | Description |
| ----------- | ----------- |
| CellDive_ID | Unique sample identifier contained in names of Halo files (images, CSV and XML boundary files)
| Patient_ID| Patient_ID from [StudyName]_Patients.xlsx corresponding to the patient from which the sample came from |
| Microscope | name of microscope used | 
| DAPI_first | string in form S[001], S[011], where the number in brackets is the number of the first staining cycle, to which all other cycles will be compared to detect cell drift/loss |
| DAPI_last | string in form S[001], S[010], where the number in brackets is the number of the last staining cycle |
Other possible columns include: Lesion_response, Treatment, Treated, etc.

<br>

#### [StudyName]_FOVs.xlsx
All meta data corresponding to individual FOVs

| Column | Description |
| ----------- | ----------- |
| CellDive_ID | Unique sample identifier of sample that FOV came from |
| FOV_number | integer representing a single FOV within a single sample (i.e., this number is unique within a sample but may be repeated across samples) |
| FOV_size | string in the form '[width]x[height]' |
| FOV_used_for_thresholding | an 'X' in this column indicates the FOV was used for manual thresholding of the sample from which it came |
| FOV_exclusion | an 'X' here indicates the FOV should NOT be included in ANY analyses |
| Marker_exclusion | comma-separated string of marker names to be excluded from an FOV |
| Num_manual_exclusions | number of manually drawn regions of cells to be excluded including epidermis, exclusion and glass regions (for QA & debugging marked/filtered exclusions) |

<br>

#### [StudyName]_Markers.xlsx
names and descriptions of all markers used in study
| Column | Description |
| ----------- | ----------- |
| Marker_name | name of marker exactly as it appears in  |
| Identity | 'Y' or 'N' indicating whether the marker is an identity marker |
| Functional | 'Y' or 'N' indicating whether the marker is a functional marker |
| Threshold_compartment | compartment used for thresholding (?) |
| Alternative_marker_name | other strings we might come across in the data that represent the same marker |
| Tumor_UMAP | 'Y' or 'N' indicating whether to use this marker in tumor UMAP |
| Immune_UMAP | 'Y' or 'N' indicating whether to use this marker in the immune UMAP |

<br>

#### [StudyName]_CellTypes.xlsx [MAJOR UPDATES 2022-10-08]
cell type definitions and labels
| Column | Description |
| --------- | -------- |
| Category | cell definition category (e.g., 'Immune', 'Tumor', 'Other', 'Functional') |
| Cell_type_full_name | full name of main cell type (e.g., 'T cell', 'Macrophage/monocyte') |
| Cell_type | abbreviated name of main cell type (e.g., 'T_All', 'Tumor_All', 'Macrophage_All') |
| Cell_type_figures | label of cell type written as it should appear in figures (may contain subscripts) |
| Cell_type_figures_subscript | if Cell_type_figures does included subscripted text, this column contains just that text (e.g., "All") |
| Subtype_full_name | full name of cell subtype (no abbreviations) |
| Subtype | abbreviated name of cell subtype (e.g., T_CD4, T_NULL, Macrophage_MHCIIpos) |
| Subtype_figures | label of cell type written as it should appear in figures (may contain formatted text such as subscripts) |
| Subtype_figures_subscript | if Subtype_figures does included subscripted text, this column contains just that text (e.g., "All") |
| Marker1 | positive or negative requirement of cell identity marker, 'Marker1' in each cell type definition (see below for details) | 
| Marker2 | positive or negative requirement of cell identity marker, 'Marker2' in each cell type definition (see below for details) |
| Marker3 | positive or negative requirement of cell identity marker, 'Marker3' in each cell type definition (see below for details) |
| MarkerN | positive or negative requirement of cell identity marker, 'MarkerN' in each cell type definition (see below for details) |

Marker columns are to be filled in as follows for basic cell type definitions:
```
no value = marker is not part of cell type definition
1 = marker must be positive to fit cell type definition
0 = marker must be negative to fit cell type definition *
```
Special cases:
It is possible to specify groups of possible positives. In some cases a cell type may be defined by one marker and/or another, or two out of three markers, for example. This is achieved by assigning a set of markers to a group identifier and appending the minimum number of positive markers to satisfy the cell type definition, in the form [GroupID_MinPositive] Examples:

| Cell type | Marker1 | Marker2 | Marker3 | Marker4 | Explanation |
| --------- | -------- | --------- | -------- | --------- | -------- |
| Cell type A | 1 | 0 | A_1 | A_1 | Marker1 and (Marker3 AND/OR Marker4) must be positive; Marker2 must be negative |
| Cell type B | A_1 | A_1 | B_1 | B_1 | (Marker1 AND/OR Marker2) and (Marker3 AND/OR Marker4) must be positive |
| Cell type C | 1 | | 0 | | Marker1 must be positive and Marker3 must be negative |


* Note: configuration may allow for flexible negative marker requirements in the case of missing markers. See script documentation for details.


<br>

#### [StudyName]_CellStates.xlsx
matrix of cell types/subtypes and the functional markers and marker combinations that should be paired to create all cell states for analysis
| Column | Description |
| --------- | -------- |
| state | functional marker name or comma-delimited string of multiple marker names |
All other columns will be either a Cell_type or Subtype from [StudyName]_CellTypes.xlsx. Fields will be marked with 'X' if the type/state pair should be made when generating cell states.

<br>

#### [StudyName]_Questions.xlsx
a systematic description of how to filter data to form sample groups for comparison
| Column | Description |
| --------- | -------- |
| Question | Biological question to be answered by this sample comparison (will be repeated two lines in a row, one line for each sample group) |
| QuestionNumber | arbitrary unique identifier of the question (will be repeated two lines in a row, one line for each sample group) |
| ComparisonVariable | the variable being examined by the question (e.g., if the question asks 'what is the difference between EBV status X and EBV status Y?' this value will be the column name from [StudyName]_Samples.xlsx or [StudyName]_FOVs.xlsx that contains EBV status.) This value will also be repeated for both groups. |
| Group | '1' or '2', indicating which group is which. Comparison will always be '2 vs 1'. |
All other column will either match a meta data variable in [StudyName]_Samples.xlsx or [StudyName]_FOVs.xlsx, 'Cell Region', 'sample_[marker]' or 'cell_[marker]'. See [ somewhere elese ] for descriptions of these variables and how to use them.

