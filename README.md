
# FastSegmentation

<!-- badges: start -->
<!-- badges: end -->

The FastSegmentation package provides a method for automatic and unsupervised cell image segmentation. 

## Installation

You can install the development version of `FastSegmentation` like so:

``` r
library(devtools)
devtools::install_github("blytheking/FastSegmentation")
```

## Segmentation Example

Segmentation by fast GP can be performed using the `generate_GP_Masks` function, with the output including the smoothed predictive mean, the binary matrix after initial, data-driven thresholding, and the final detected cell masks.

``` r
## Perform segmentation
library(FastSegmentation)

#Set path to image to segment
file_path <- "path/to/file/img.jpg"

#Run segmentation function
gp_masks_result <- generate_GP_Masks(file_path)

## Visualize results
library(plot3D)
library(magick)

#Original image
img <- image_read(file_path)
ori_img_matrix <- as.numeric(img[[1]])[,,1]
image2D(ori_img_matrix, main = "Original Image")

#Predictive Mean
image2D(gp_masks_result$combined_predmean, main = "Predictive Mean")

#Thresholding by Criterion 1
image2D(gp_masks_result$combined_thresholded1, main = "Binary Matrix")

#Final Cell Masks
image2D(gp_masks_result, main = "Binary Matrix")
```

## Notes
`FastSegmentation` uses the `EBImage` package from `BiocManager`. Since `EBImage` is not part of `CRAN`, please first install it using this code:

``` r
library(BiocManager)
BiocManager::install("EBImage")
library(EBImage)
```

## Relevant Literature
This package was adapted from the methods described in the following manuscript:

Baracaldo, L., King, B., Yan, H., Lin, Y., Miolane, N., & Gu, M. (2025). Unsupervised cell segmentation by fast Gaussian processes. _arXiv preprint arXiv:2505.18902_.

