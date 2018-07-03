RUV-III Challenges
================

-   [Challenge 1 data](#challenge-1-data)
-   [Challenge 2 data](#challenge-2-data)

The two challenges below are based on the RUV method for normalisation, RUV-III. There are two datasets provided containing gene expression data: `Challenge1Data.csv` and `Challenge2Data.csv`. The dataset `Challenge.information.csv` contains clinical information of the samples.

There are 5 different biological condition of interest (`Challenge.information.csv` column `BiologicalCondition`) which have been run in 5 different batches.

The two datasets are the same but one has already been normalised and the other one hasn't. **Tasks**

Use `Challenge1Data.csv` or `Challenge2Data.csv` and assess if:

-   Data needs normalisation
-   If the data needs normalisation, explain why and apply a normalisation.
-   If you think the data has already been normalised, explain why and comment on the quality of the normalisation.
-   The same set of negative control genes is provided in both challenges. It is up to you to decide whether you need it!

**Hints**

-   Remember that there is not one way to tell if the data needs normalisation
-   Remember to check the biology! Consider analysing the expression of the Y chromosome genes `ZFY` and `EIF1AY` to make your decision. These genes on the Y chromosome should be expressed in male samples but not in female samples.

Challenge 1 data
----------------

``` r
ch1_data <- read.csv("Challenge1Data.csv", row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)

Con100 <- readRDS('controlGene100.rds')

Info <- read.csv('Challenge.information.csv', row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)

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

| SampleIds | BiologicalConditions | Gender |
|:----------|:---------------------|:-------|
| Sample\_1 | Conditions A         | female |
| Sample\_2 | Conditions A         | female |
| Sample\_3 | Conditions A         | female |
| Sample\_4 | Conditions A         | female |
| Sample\_5 | Conditions A         | female |
| Sample\_6 | Conditions A         | male   |

Challenge 2 data
----------------

``` r
ch2_data <- read.csv("Challenge2Data.csv", row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)
Con100 <- readRDS('controlGene100.rds')
Info <- read.csv('Challenge.information.csv', row.names = 1, as.is = TRUE, stringsAsFactors = FALSE)

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

| SampleIds | BiologicalConditions | Gender |
|:----------|:---------------------|:-------|
| Sample\_1 | Conditions A         | female |
| Sample\_2 | Conditions A         | female |
| Sample\_3 | Conditions A         | female |
| Sample\_4 | Conditions A         | female |
| Sample\_5 | Conditions A         | female |
| Sample\_6 | Conditions A         | male   |
