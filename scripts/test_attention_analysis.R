# TODO: this does not work. need to revise.


#' Test Suite for Attention Analysis Functions (Base R version)
#'
#' Self-contained tests using synthetic data to verify the attention
#' heatmap pipeline works correctly end-to-end.
#'
#' Run with: Rscript test-attention-base.R

# Source the attention analysis functions
# source("attention-analysis-base.R")
devtools::load_all()

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("ATTENTION ANALYSIS TEST SUITE (Base R)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# =============================================================================
# Synthetic Data Generators
# =============================================================================

generate_synthetic_gaze <- function(n_fixations = 10,
                                    samples_per_fixation = 30,
                                    saccade_samples = 5,
                                    sampling_rate = 60,
                                    fixation_noise_sd = 0.01) {

  set.seed(42)

  # Generate fixation centers
  fixation_centers <- data.frame(
    fix_id = 1:n_fixations,
    center_x = runif(n_fixations, 0.2, 0.8),
    center_y = runif(n_fixations, 0.2, 0.8)
  )

  # Build gaze stream
  gaze_list <- list()
  current_time <- 0
  sample_id <- 1

  for (i in 1:n_fixations) {
    # Fixation period
    fixation_data <- data.frame(
      sample_id = sample_id:(sample_id + samples_per_fixation - 1),
      time_sec = current_time + (0:(samples_per_fixation - 1)) / sampling_rate,
      gaze_x = fixation_centers$center_x[i] + rnorm(samples_per_fixation, 0, fixation_noise_sd),
      gaze_y = fixation_centers$center_y[i] + rnorm(samples_per_fixation, 0, fixation_noise_sd),
      true_fixation_id = i,
      is_true_fixation = TRUE
    )

    gaze_list[[length(gaze_list) + 1]] <- fixation_data
    sample_id <- sample_id + samples_per_fixation
    current_time <- current_time + samples_per_fixation / sampling_rate

    # Saccade to next fixation
    if (i < n_fixations) {
      saccade_x <- seq(fixation_centers$center_x[i],
                       fixation_centers$center_x[i + 1],
                       length.out = saccade_samples)
      saccade_y <- seq(fixation_centers$center_y[i],
                       fixation_centers$center_y[i + 1],
                       length.out = saccade_samples)

      saccade_data <- data.frame(
        sample_id = sample_id:(sample_id + saccade_samples - 1),
        time_sec = current_time + (0:(saccade_samples - 1)) / sampling_rate,
        gaze_x = saccade_x,
        gaze_y = saccade_y,
        true_fixation_id = NA,
        is_true_fixation = FALSE
      )

      gaze_list[[length(gaze_list) + 1]] <- saccade_data
      sample_id <- sample_id + saccade_samples
      current_time <- current_time + saccade_samples / sampling_rate
    }
  }

  gaze_data <- do.call(rbind, gaze_list)
  gaze_data$plane <- "AXIAL"
  gaze_data$slice_index <- 0.5
  gaze_data$time_aligned <- gaze_data$time_sec

  attr(gaze_data, "n_true_fixations") <- n_fixations
  return(gaze_data)
}

generate_mock_nifti <- function(dims = c(64, 64, 20)) {
  vol <- array(runif(prod(dims)) * 255, dim = dims)

  mock_nvimage <- list(
    matRAS = diag(4),
    dimsRAS = c(3, dims),
    pixDimsRAS = c(1, 1, 1, 5)
  )

  return(list(
    data = vol,
    dims = dims,
    nvimage = mock_nvimage,
    nifti = NULL
  ))
}

# Mock function for testing (filters by plane)
get_all_gaze_for_slice <- function(integrated_data, slice_num, plane_filter, dims) {
  integrated_data[integrated_data$plane == plane_filter, ]
}

# =============================================================================
# Run Tests
# =============================================================================

passed <- 0
failed <- 0

# Test 1: Fixation Detection
cat("Test 1: Fixation Detection\n")
tryCatch({
  gaze_data <- generate_synthetic_gaze(n_fixations = 10, fixation_noise_sd = 0.002)
  # Use higher threshold - synthetic data at 60Hz with 0.002 noise SD gives ~0.17 velocity
  result <- detect_fixations(gaze_data, velocity_threshold = 0.5, min_fixation_duration = 0.05)
  n_detected <- length(unique(na.omit(result$fixation_id)))

  if (n_detected >= 8 && n_detected <= 12) {
    cat("  ✓ PASS: Detected", n_detected, "fixations (expected ~10)\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Detected", n_detected, "fixations (expected ~10)\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 2: Attention Weights
cat("\nTest 2: Attention Weights\n")
tryCatch({
  gaze_data <- generate_synthetic_gaze(n_fixations = 5, fixation_noise_sd = 0.002)
  with_fixations <- detect_fixations(gaze_data, velocity_threshold = 0.5, min_fixation_duration = 0.05)
  with_weights <- calculate_attention_weights(with_fixations)

  fix_wt <- mean(with_weights$attention_weight[with_weights$is_fixation], na.rm = TRUE)
  sac_wt <- mean(with_weights$attention_weight[!with_weights$is_fixation], na.rm = TRUE)

  if (fix_wt > sac_wt && all(with_weights$attention_weight >= 0, na.rm = TRUE)) {
    cat("  ✓ PASS: Fixation weight (", round(fix_wt, 3), ") > saccade weight (", round(sac_wt, 3), ")\n", sep="")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Weight comparison failed\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 3: Weighted KDE Output Dimensions
cat("\nTest 3: Weighted KDE Output Dimensions\n")
tryCatch({
  x <- c(0.3, 0.3, 0.7, 0.7)
  y <- c(0.3, 0.7, 0.3, 0.7)
  weights <- c(1, 1, 1, 1)

  result <- weighted_kde2d(x, y, weights, n = c(32, 32), bandwidth = 0.1, xlim = c(0, 1), ylim = c(0, 1))

  if (all(dim(result$z) == c(32, 32)) && max(result$z) > 0) {
    cat("  ✓ PASS: Output dimensions 32x32, max density =", round(max(result$z), 4), "\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Dimension or density issue\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 4: Weight Sensitivity in KDE
cat("\nTest 4: Weight Sensitivity in KDE\n")
tryCatch({
  x <- c(rep(0.25, 10), rep(0.75, 10))
  y <- c(rep(0.5, 10), rep(0.5, 10))
  weights <- c(rep(10, 10), rep(1, 10))

  result <- weighted_kde2d(x, y, weights, n = c(50, 50), bandwidth = 0.05, xlim = c(0, 1), ylim = c(0, 1))

  left_idx <- which.min(abs(result$x - 0.25))
  right_idx <- which.min(abs(result$x - 0.75))
  center_y_idx <- which.min(abs(result$y - 0.5))

  density_left <- result$z[left_idx, center_y_idx]
  density_right <- result$z[right_idx, center_y_idx]

  if (density_left > density_right) {
    cat("  ✓ PASS: High-weight cluster (", round(density_left, 4), ") > low-weight (", round(density_right, 4), ")\n", sep="")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Weights not respected\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 5: Slice Heatmap Generation
cat("\nTest 5: Slice Heatmap Generation\n")
tryCatch({
  gaze_data <- generate_synthetic_gaze(n_fixations = 8, fixation_noise_sd = 0.002)
  nifti_data <- generate_mock_nifti(dims = c(64, 64, 20))

  result <- generate_slice_heatmap(
    integrated_data = gaze_data,
    nifti_data = nifti_data,
    slice_num = 10,
    plane = "AXIAL",
    bandwidth = 3,
    normalize = TRUE
  )

  if (is.matrix(result$heatmap) && max(result$heatmap) <= 1 && max(result$heatmap) > 0) {
    cat("  ✓ PASS: Heatmap", nrow(result$heatmap), "x", ncol(result$heatmap),
        ", n_samples =", result$n_samples, ", max =", round(max(result$heatmap), 3), "\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Heatmap generation issue\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 6: Volume Heatmap Generation
cat("\nTest 6: Volume Heatmap Generation\n")
tryCatch({
  gaze_data <- generate_synthetic_gaze(n_fixations = 5, fixation_noise_sd = 0.002)
  nifti_data <- generate_mock_nifti(dims = c(32, 32, 10))

  result <- generate_volume_heatmap(
    integrated_data = gaze_data,
    nifti_data = nifti_data,
    plane = "AXIAL",
    bandwidth = 2,
    verbose = FALSE
  )

  if (all(dim(result$volume) == c(32, 32, 10)) && max(result$volume) <= 1) {
    cat("  ✓ PASS: Volume", paste(dim(result$volume), collapse="x"),
        ", slices with data =", result$n_slices_with_data, "\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Volume generation issue\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 7: Spatial Accuracy
cat("\nTest 7: Spatial Accuracy - Attention Concentrates Correctly\n")
tryCatch({
  set.seed(123)
  # Create fixation-like data: clustered with small movements
  n_samples <- 200
  gaze_x <- numeric(n_samples)
  gaze_y <- numeric(n_samples)

  # 5 fixation clusters in upper-right region
  centers_x <- c(0.65, 0.70, 0.75, 0.80, 0.85)
  centers_y <- c(0.65, 0.70, 0.75, 0.80, 0.85)
  samples_per_fix <- 40

  for (i in 1:5) {
    idx <- ((i-1)*samples_per_fix + 1):(i*samples_per_fix)
    gaze_x[idx] <- centers_x[i] + rnorm(samples_per_fix, 0, 0.005)
    gaze_y[idx] <- centers_y[i] + rnorm(samples_per_fix, 0, 0.005)
  }

  gaze_data <- data.frame(
    sample_id = 1:n_samples,
    time_sec = (1:n_samples) / 60,
    gaze_x = gaze_x,
    gaze_y = gaze_y,
    plane = "AXIAL",
    slice_index = 0.5,
    time_aligned = (1:n_samples) / 60
  )
  nifti_data <- generate_mock_nifti(dims = c(64, 64, 20))

  result <- generate_slice_heatmap(gaze_data, nifti_data, slice_num = 10, plane = "AXIAL", bandwidth = 3)

  heatmap <- result$heatmap
  upper_right <- mean(heatmap[33:64, 33:64])
  lower_left <- mean(heatmap[1:32, 1:32])
  ratio <- upper_right / max(lower_left, 1e-10)

  if (upper_right > lower_left * 2) {
    cat("  ✓ PASS: Attention concentrated in correct quadrant (ratio =", round(ratio, 1), ")\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Attention not concentrated in correct region\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# Test 8: CSV Export
cat("\nTest 8: CSV Export\n")
tryCatch({
  gaze_data <- generate_synthetic_gaze(n_fixations = 5, fixation_noise_sd = 0.002)
  nifti_data <- generate_mock_nifti(dims = c(32, 32, 10))

  volume_result <- generate_volume_heatmap(gaze_data, nifti_data, plane = "AXIAL", bandwidth = 2, verbose = FALSE)

  temp_file <- tempfile(fileext = ".csv")
  export_attention_csv(volume_result, temp_file, threshold = 0.01)

  exported <- read.csv(temp_file)
  unlink(temp_file)

  if (nrow(exported) > 0 && all(c("voxel_x", "voxel_y", "voxel_z", "attention") %in% names(exported))) {
    cat("  ✓ PASS: Exported", nrow(exported), "voxels, top attention =", round(max(exported$attention), 3), "\n")
    passed <- passed + 1
  } else {
    cat("  ✗ FAIL: Export issue\n")
    failed <- failed + 1
  }
}, error = function(e) {
  cat("  ✗ ERROR:", e$message, "\n")
  failed <<- failed + 1
})

# =============================================================================
# Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("RESULTS:", passed, "passed,", failed, "failed\n")
if (failed == 0) {
  cat("✅ ALL TESTS PASSED!\n")
} else {
  cat("❌ Some tests failed\n")
}
cat(paste(rep("=", 60), collapse = ""), "\n")
