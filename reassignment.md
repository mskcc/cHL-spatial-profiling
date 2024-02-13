# Hodgkin cell reassignment hierarchy

*completed at the FOV-level*

*steps 3 and 4 are interchangeable*

Step 1: Recovery of HRS cells

1. Take UNKNOWN population
2. Drop CD30 completeness threshold for cells that meet CD30 intensity threshold
3. Drop MUM1 completeness threshold for cells that meet MUM1 intensity threshold
    1. Some cells become an HRS cell

*FOVs with <0.1% HRS cells (over DAPI) following Step 1 are removed from the cohort*

Step 2: Resolution of CD3/CD20 double positive cells

1. Take UNKNOWN population
2. Rank reassignment for CD3/CD20 (0% CD3/CD20 double positive allowed)
    1. If CD20 “wins”: CD3 positivity switches to negative
        1. Any CD4, CD8, FOXP3 positivity switches to negative
        2. Still potential for some UNKNOWN cells based on positivity of other markers       
    2. If CD3 “wins”: CD20 positivity switches to negative
        1. Some cells will have new T cell subtype identity
        2. Still potential for some UNKNOWN cells based on positivity of other markers


Step 3-A: Resolution of T cell/null phenotype population

0.  *MUST BE DONE AFTER Step 2 and after type reassignment*; *5% T cell/null phenotype cells over total T cell population allowed*
1.  Take T cell/null phenotype population (double neg cells CD4/CD8)
2.  If PCT[CD4-,CD8-]>5% then drop CD4,CD8 completeness threshold for cells that meet cytoplasmic intensity threshold
3.  If PCT[CD4-,CD8-]>5% then take remainder of cells (cells that both fail the CD4 and CD8 intensity
    thresholds)
    1.  Rank reassignment for CD4/CD8 (5% T cell/null phenotype cells over total T cell population allowed)
        1. If CD4 “wins”: CD4 positivity switches to positive
            - Cell becomes a CD4+ T cell
        1. If CD8 “wins”: CD8 positivity switches to positive
            - Cell becomes a CD8+ T cell

Step 3-B: Resolution of CD4+/CD8+ T cells

0.  *MUST BE DONE AFTER Step 2, Step 3-A and after type reassignment*; *5% T CD4+/CD8+ phenotype cells over total T cell population allowed*
1.  Take CD4+/CD8+ T cell population
2.  Rank assignment for CD4/CD8 (5% CD4+/CD8+ T cells over total T cell population allowed)
    1.  If CD4 “wins”: CD8 positivity switches to negative
        1. Cell becomes a CD4+ T cell or CD4+ regulatory T cell
    2.  If CD8 “wins”: CD4 positivity switches to negative
        1. Cell becomes a CD8+ T cell or CD8+ regulatory T cell


