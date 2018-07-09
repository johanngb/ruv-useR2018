RUV-III Challenges
================

-   [Tasks](#tasks)
-   [Challenge 1 data](#challenge-1-data)
-   [Challenge 2 data](#challenge-2-data)

The two challenges below are based on the RUV method for normalisation, RUV-III. There are two datasets provided containing gene expression data obtained with Nanostring technology: `Challenge1Data.csv` and `Challenge2Data.csv`. The dataset `Challenge.information.csv` contains the clinical information of the samples (biological conditions of interest and gender of the samples).

There are 5 different biological conditions of interest (in `Challenge.information.csv` column `Biology`) which involve T-cells treated in different ways. They have been run in 2 different runs: `Run1` and `Run2`.

Tasks
-----

The two datasets are the same but one has already been normalised and the other one hasn't. Use `Challenge1Data.csv` or `Challenge2Data.csv` and assess if:

-   The gene expression needs normalisation
-   If the data needs normalisation, explain why and apply a normalisation.
-   If you think the data has already been normalised, explain how you deduce it and comment on the quality of the normalisation.
-   The same set of negative control genes is provided in both challenges. It is up to you to decide whether you need it!

**Hints**

-   Remember that there is not only one way to tell if the data needs normalisation
-   Remember to check the known biology! Consider analysing the expression of the Y chromosome genes `ZFY` and `EIF1AY` to inform your decision. These genes on the Y chromosome should be expressed in male samples but not in female samples.

Challenge 1 data
----------------

``` r
ch1_data <- read.csv("Challenge1Data.csv", row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)

Con100 <- readRDS('controlGene100.rds')

Info <- read.csv('Challenge.information.csv', as.is = TRUE, stringsAsFactors = FALSE)


ch1_data <- as.matrix(ch1_data)
kable(ch1_data[1:4,1:4])
```

|         |  Sample\_1|  Sample\_2|  Sample\_3|  Sample\_4|
|---------|----------:|----------:|----------:|----------:|
| ACTB    |  12.051549|  11.912141|  11.956739|  11.349281|
| AHR     |  10.049848|   9.902375|   9.961450|   9.290019|
| AIF1    |   3.000000|   2.000000|   3.700440|   3.321928|
| ALDH7A1 |   1.584963|   1.000000|   3.169925|   2.807355|

``` r
kable(head(Info))
```

| SampleIds | Gender | Runs | Biology           |
|:----------|:-------|:-----|:------------------|
| Sample\_1 | female | Run1 | Activated\_4hours |
| Sample\_2 | female | Run1 | Activated\_4hours |
| Sample\_3 | female | Run1 | Activated\_4hours |
| Sample\_4 | female | Run1 | Activated\_4hours |
| Sample\_5 | female | Run1 | Activated\_4hours |
| Sample\_6 | male   | Run1 | Activated\_4hours |

Challenge 2 data
----------------

``` r
ch2_data <- read.csv("Challenge2Data.csv", row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)

Con100 <- readRDS('controlGene100.rds')
Info <- read.csv('Challenge.information.csv', as.is = TRUE, stringsAsFactors = FALSE)

ch2_data <- as.matrix(ch2_data)
kable(ch2_data[1:4,1:4])
```

|         |   Sample\_1|   Sample\_2|  Sample\_3|  Sample\_4|
|---------|-----------:|-----------:|----------:|----------:|
| ACTB    |   4.0638266|   3.6081448|  3.6605152|  3.8799500|
| AHR     |   2.4416056|   1.8788088|  2.0303678|  1.9730563|
| AIF1    |   0.9261947|   0.2359753|  1.4106417|  2.1742318|
| ALDH7A1 |  -0.5912845|  -0.8199275|  0.8812682|  0.9607422|

``` r
kable(head(Info))
```

| SampleIds | Gender | Runs | Biology           |
|:----------|:-------|:-----|:------------------|
| Sample\_1 | female | Run1 | Activated\_4hours |
| Sample\_2 | female | Run1 | Activated\_4hours |
| Sample\_3 | female | Run1 | Activated\_4hours |
| Sample\_4 | female | Run1 | Activated\_4hours |
| Sample\_5 | female | Run1 | Activated\_4hours |
| Sample\_6 | male   | Run1 | Activated\_4hours |
