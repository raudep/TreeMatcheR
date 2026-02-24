
# TreeMatcheR

`TreeMatcheR` provides tools for the co-registration of single tree data
in forest stands and forest plots. It is applicable to both static and
dynamic data capture (e.g., SLAM-based moving sensors).

The package implements a stem diameter weighted linking algorithm that
improves linking accuracy when operating on diverse diameter stands with
stem position errors. This package implements the methodology described
in Olofsson and Holmgren (2022) <doi:10.14214/sf.10712>.

## Installation

You can install the development version of TreeMatcheR from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("raudep/TreeMatcheR")
```

## Example Usage

Here is a basic example showing how to run the built-in demonstration.
The algorithm requires defining column headers so it knows where to find
coordinates (X, Y) and stem size (R) for the weighting process.

``` r
library(TreeMatcheR)

# Run the built-in demo process
runDemoPlotMatch()
```

### Manual Data Organization

To run this on your own dataframes, organize your tables and map your
column names using the `TableHeaderLabels()` function.

``` r
# 1. Map your specific column names
metaAir <- TableHeaderLabels(
  id = "ID", x = "X", y = "Y", z = NULL, 
  r = "H", t = NULL, scaleR = 0.0075, scan = NULL
)

metaGround <- TableHeaderLabels(
  id = c("PLOT", "ID"), x = "X", y = "Y", z = NULL, 
  r = "STEMDIAM", t = NULL, scaleR = 0.005, scan = NULL
)

# 2. Run the high-performance matching algorithm
# searchDist: Radius in meters to look for candidate matches. 
# estX0/estY0: Field plot center coordinates (should be in the same coordinate system as the the airborne data - air_data)
result <- tableMatchAndGetLinkTableAndTrafoMatLocal(
  tbForest = air_data, 
  tbPlot = ground_data, 
  metaForest = metaAir, 
  metaPlot = metaGround,
  searchDist = 5.0, 
  estX0 = 311.99, 
  estY0 = 212.44
)

# 3. Access the resulting transformation matrix
print(result$M)
```
