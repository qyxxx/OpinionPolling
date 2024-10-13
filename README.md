# A Model-free Framework for Evaluating the Reliability of a New Device with Multiple Imperfect Reference Standards

A common practice for establishing the reliability of a new computer-aided diagnostic (CAD) device is to evaluate how well its clinical measurements agree with those of a gold standard test. However, in many clinical studies, a gold standard is unavailable, and one needs to aggregate information from multiple imperfect reference standards for evaluation. A key challenge here is the heterogeneity in diagnostic accuracy across different reference standards, which may lead to biased evaluation of a device if improperly accounted for during the aggregation process. We propose an intuitive and easy-to-use statistical framework for evaluation of a device by assessing agreement between its measurements and the weighted sum of measurements from multiple imperfect reference standards, where weights representing relative reliability of each reference standard are determined by a model-free, unsupervised inductive procedure. Specifically, the inductive procedure recursively assigns higher weights to reference standards whose assessments are more consistent with each other and form a majority opinion, while assigning lower weights to those with greater discrepancies. Unlike existing methods, our approach does not require any modeling assumptions or external data to quantify heterogeneous accuracy levels of reference standards. It only requires specifying an appropriate agreement index used for weight assignment and device evaluation. The framework is applied to evaluate a CAD device for kidney obstruction by comparing its diagnostic ratings with those of multiple nuclear medicine physicians.

## Installation

You can install the development version of SLCARE like so:

``` r
if (!require("pak", quietly = TRUE))
    install.packages("pak")

pak::pak("qyxxx/OpinionPolling")
```
