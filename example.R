##########################################################################
## example.R
## 
## FastSegmentation Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2026 - present Blythe King
##							  
##    
##########################################################################

################################################################
# Example segmentation code for microscopy images
################################################################

library(FastSegmentation)
library(plot3D)
library(magick)

#------------------------------------------------------------------------------
# Example 1: Full Segmentation Workflow for Nuclear Channel
#------------------------------------------------------------------------------

# File path to TIF, PNG, or JPEG image
img_path <- system.file("extdata", "example_cells.jpg", package = "FastSegmentation")

# ---- Run segmentation ----
gp_masks_result <- generate_GP_Masks(img_path)

# ---- Visualization ----

# Original image
img <- image_read(img_path)
ori_img_matrix <- as.numeric(img[[1]])[,,1]
image2D(ori_img_matrix, main = "Original Image", axes=FALSE, xlab="", ylab="")

# Predictive Mean
image2D(gp_masks_result$combined_predmean, main = "Predictive Mean", 
        axes=FALSE, xlab="", ylab="")

# Thresholding by Criterion 1
image2D(gp_masks_result$combined_thresholded1, 
        main = "Thresholded Image", axes=FALSE, xlab="", ylab="")

# Final Cell Masks
image2D(gp_masks_result$GP_masks, main = "Cell Masks", axes=FALSE, xlab="", ylab="")

#------------------------------------------------------------------------------
# Example 2: Full Segmentation Workflow for Whole Cell Channel
#------------------------------------------------------------------------------

# File path to TIF, PNG, or JPEG image
img_path <- system.file("extdata", "whole_cell_example.tif",
                        package = "FastSegmentation")

# ---- Run segmentation ----
gp_masks_result <- generate_GP_Masks(img_path)

# ---- Visualize ----

# Original image
img <- image_read(img_path)
ori_img_matrix <- as.numeric(img[[1]])[,,1]
image2D(ori_img_matrix, main = "Original Image", axes=FALSE, xlab="", ylab="")

# Predictive Mean
image2D(gp_masks_result$combined_predmean, main = "Predictive Mean", 
        axes=FALSE, xlab="", ylab="")

# Thresholding by Criterion 1
image2D(gp_masks_result$combined_thresholded1, 
        main = "Thresholded Image", axes=FALSE, xlab="", ylab="")

# Final Cell Masks
image2D(gp_masks_result$GP_masks, main = "Cell Masks", axes=FALSE, xlab="", ylab="")

