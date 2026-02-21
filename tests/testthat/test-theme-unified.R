# Unit Tests for Unified Theme System
# Tests for theme_unified.R functions

library(testthat)

# Skip all tests if ggplot2 is not available
skip_if_not_installed("ggplot2")
library(ggplot2)

# Test theme_mixsynergy function
test_that("theme_mixsynergy returns a valid ggplot2 theme", {
  theme_obj <- theme_mixsynergy()
  
  # Check that it returns a theme object
  expect_s3_class(theme_obj, "theme")
  expect_s3_class(theme_obj, "gg")
  
  # Check that it has expected components
  expect_true(!is.null(theme_obj$plot.title))
  expect_true(!is.null(theme_obj$axis.title))
  expect_true(!is.null(theme_obj$legend.title))
})

test_that("theme_mixsynergy accepts base_size parameter", {
  theme_small <- theme_mixsynergy(base_size = 8)
  theme_large <- theme_mixsynergy(base_size = 16)
  
  # Both should be valid themes
  expect_s3_class(theme_small, "theme")
  expect_s3_class(theme_large, "theme")
  
  # Themes with different base sizes should be different
  expect_false(identical(theme_small, theme_large))
})

test_that("theme_mixsynergy accepts base_family parameter", {
  theme_default <- theme_mixsynergy(base_family = "")
  theme_custom <- theme_mixsynergy(base_family = "sans")
  
  # Both should be valid themes
  expect_s3_class(theme_default, "theme")
  expect_s3_class(theme_custom, "theme")
})

test_that("theme_mixsynergy can be applied to a ggplot", {
  # Create a simple plot
  df <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot(df, aes(x = x, y = y)) + 
    geom_point() + 
    theme_mixsynergy()
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the theme was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test theme_mixsynergy_heatmap function
test_that("theme_mixsynergy_heatmap returns a valid ggplot2 theme", {
  theme_obj <- theme_mixsynergy_heatmap()
  
  # Check that it returns a theme object
  expect_s3_class(theme_obj, "theme")
  expect_s3_class(theme_obj, "gg")
  
  # Check that it has expected components
  expect_true(!is.null(theme_obj$panel.grid.major))
  expect_true(!is.null(theme_obj$panel.grid.minor))
})

test_that("theme_mixsynergy_heatmap accepts base_size parameter", {
  theme_small <- theme_mixsynergy_heatmap(base_size = 8)
  theme_large <- theme_mixsynergy_heatmap(base_size = 16)
  
  # Both should be valid themes
  expect_s3_class(theme_small, "theme")
  expect_s3_class(theme_large, "theme")
})

test_that("theme_mixsynergy_heatmap can be applied to a heatmap", {
  # Create a simple heatmap
  df <- expand.grid(x = 1:5, y = 1:5)
  df$value <- rnorm(25)
  
  p <- ggplot(df, aes(x = x, y = y, fill = value)) + 
    geom_tile() + 
    theme_mixsynergy_heatmap()
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the theme was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test scale_fill_mixsynergy function
test_that("scale_fill_mixsynergy returns a valid scale", {
  scale_obj <- scale_fill_mixsynergy()
  
  # Check that it returns a scale object
  expect_s3_class(scale_obj, "Scale")
  expect_s3_class(scale_obj, "ScaleContinuous")
})

test_that("scale_fill_mixsynergy accepts option parameter", {
  scale_viridis <- scale_fill_mixsynergy(option = "viridis")
  scale_magma <- scale_fill_mixsynergy(option = "magma")
  scale_plasma <- scale_fill_mixsynergy(option = "plasma")
  
  # All should be valid scales
  expect_s3_class(scale_viridis, "Scale")
  expect_s3_class(scale_magma, "Scale")
  expect_s3_class(scale_plasma, "Scale")
})

test_that("scale_fill_mixsynergy accepts direction parameter", {
  scale_forward <- scale_fill_mixsynergy(direction = 1)
  scale_reverse <- scale_fill_mixsynergy(direction = -1)
  
  # Both should be valid scales
  expect_s3_class(scale_forward, "Scale")
  expect_s3_class(scale_reverse, "Scale")
})

test_that("scale_fill_mixsynergy can be applied to a plot", {
  # Create a simple plot with continuous fill
  df <- expand.grid(x = 1:5, y = 1:5)
  df$value <- rnorm(25)
  
  p <- ggplot(df, aes(x = x, y = y, fill = value)) + 
    geom_tile() + 
    scale_fill_mixsynergy()
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the scale was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test scale_fill_mixsynergy_diverging function
test_that("scale_fill_mixsynergy_diverging returns a valid scale", {
  scale_obj <- scale_fill_mixsynergy_diverging()
  
  # Check that it returns a scale object
  expect_s3_class(scale_obj, "Scale")
  expect_s3_class(scale_obj, "ScaleContinuous")
})

test_that("scale_fill_mixsynergy_diverging accepts midpoint parameter", {
  scale_zero <- scale_fill_mixsynergy_diverging(midpoint = 0)
  scale_custom <- scale_fill_mixsynergy_diverging(midpoint = 0.5)
  
  # Both should be valid scales
  expect_s3_class(scale_zero, "Scale")
  expect_s3_class(scale_custom, "Scale")
})

test_that("scale_fill_mixsynergy_diverging accepts limits parameter", {
  scale_default <- scale_fill_mixsynergy_diverging()
  scale_custom <- scale_fill_mixsynergy_diverging(limits = c(-1, 1))
  
  # Both should be valid scales
  expect_s3_class(scale_default, "Scale")
  expect_s3_class(scale_custom, "Scale")
})

test_that("scale_fill_mixsynergy_diverging can be applied to a plot", {
  # Create a simple plot with diverging fill
  df <- expand.grid(x = 1:5, y = 1:5)
  df$value <- rnorm(25, mean = 0, sd = 1)
  
  p <- ggplot(df, aes(x = x, y = y, fill = value)) + 
    geom_tile() + 
    scale_fill_mixsynergy_diverging()
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the scale was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test scale_color_mixsynergy function
test_that("scale_color_mixsynergy returns a valid scale", {
  scale_obj <- scale_color_mixsynergy()
  
  # Check that it returns a scale object
  expect_s3_class(scale_obj, "Scale")
  expect_s3_class(scale_obj, "ScaleDiscrete")
})

test_that("scale_color_mixsynergy accepts option parameter", {
  scale_viridis <- scale_color_mixsynergy(option = "viridis")
  scale_magma <- scale_color_mixsynergy(option = "magma")
  
  # Both should be valid scales
  expect_s3_class(scale_viridis, "Scale")
  expect_s3_class(scale_magma, "Scale")
})

test_that("scale_color_mixsynergy can be applied to a plot", {
  # Create a simple plot with discrete color
  df <- data.frame(
    x = 1:20,
    y = rnorm(20),
    group = rep(c("A", "B", "C", "D"), 5)
  )
  
  p <- ggplot(df, aes(x = x, y = y, color = group)) + 
    geom_point() + 
    scale_color_mixsynergy()
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the scale was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Test facet_mixsynergy function
test_that("facet_mixsynergy returns a valid facet specification", {
  facet_obj <- facet_mixsynergy(rows = "Cluster", cols = "Factor")
  
  # Check that it returns a facet object
  expect_s3_class(facet_obj, "FacetGrid")
  expect_s3_class(facet_obj, "Facet")
})

test_that("facet_mixsynergy accepts rows parameter only", {
  facet_obj <- facet_mixsynergy(rows = "Cluster")
  
  # Check that it returns a facet object
  expect_s3_class(facet_obj, "FacetGrid")
})

test_that("facet_mixsynergy accepts cols parameter only", {
  facet_obj <- facet_mixsynergy(cols = "Factor")
  
  # Check that it returns a facet object
  expect_s3_class(facet_obj, "FacetGrid")
})

test_that("facet_mixsynergy requires at least one parameter", {
  # Should error when neither rows nor cols is specified
  expect_error(
    facet_mixsynergy(),
    "At least one of 'rows' or 'cols' must be specified"
  )
})

test_that("facet_mixsynergy accepts scales parameter", {
  facet_fixed <- facet_mixsynergy(rows = "Cluster", cols = "Factor", scales = "fixed")
  facet_free <- facet_mixsynergy(rows = "Cluster", cols = "Factor", scales = "free")
  facet_free_y <- facet_mixsynergy(rows = "Cluster", cols = "Factor", scales = "free_y")
  
  # All should be valid facet objects
  expect_s3_class(facet_fixed, "FacetGrid")
  expect_s3_class(facet_free, "FacetGrid")
  expect_s3_class(facet_free_y, "FacetGrid")
})

test_that("facet_mixsynergy can be applied to a plot", {
  # Create a simple plot with faceting
  df <- data.frame(
    x = rep(1:10, 4),
    y = rnorm(40),
    Cluster = rep(c(1, 1, 2, 2), each = 10),
    Factor = rep(c("F1", "F2"), each = 20)
  )
  
  p <- ggplot(df, aes(x = x, y = y)) + 
    geom_line() + 
    facet_mixsynergy(rows = "Cluster", cols = "Factor")
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the facet was applied
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

# Integration tests: Test that all components work together
test_that("unified theme system components work together", {
  # Create a comprehensive plot using all theme components
  df <- expand.grid(
    x = 1:5,
    y = 1:5,
    Cluster = c("C1", "C2")
  )
  df$value <- rnorm(50)
  
  p <- ggplot(df, aes(x = x, y = y, fill = value)) + 
    geom_tile() + 
    scale_fill_mixsynergy_diverging() +
    facet_mixsynergy(cols = "Cluster") +
    theme_mixsynergy_heatmap() +
    labs(title = "Test Heatmap", x = "X Axis", y = "Y Axis")
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the plot builds successfully
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

test_that("unified theme system works with line plots", {
  # Create a line plot using theme components
  df <- data.frame(
    Time = rep(1:20, 4),
    Score = rnorm(80),
    Cluster = rep(c(1, 2), each = 40),
    Factor = rep(c("F1", "F2"), each = 20)
  )
  
  p <- ggplot(df, aes(x = Time, y = Score, group = Cluster, color = factor(Cluster))) + 
    geom_line() + 
    scale_color_mixsynergy() +
    facet_mixsynergy(cols = "Factor", scales = "free_y") +
    theme_mixsynergy() +
    labs(title = "Test Line Plot", x = "Time", y = "Score")
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the plot builds successfully
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})

test_that("unified theme system works with bar plots", {
  # Create a bar plot using theme components
  df <- data.frame(
    Muscle = rep(paste0("M", 1:5), 2),
    Loading = rnorm(10),
    Cluster = rep(c("C1", "C2"), each = 5)
  )
  
  p <- ggplot(df, aes(x = Muscle, y = Loading, fill = Cluster)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_color_mixsynergy() +
    facet_mixsynergy(cols = "Cluster") +
    theme_mixsynergy() +
    labs(title = "Test Bar Plot", x = "Muscle", y = "Loading")
  
  # Check that plot can be built without errors
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
  
  # Check that the plot builds successfully
  built <- ggplot_build(p)
  expect_true(!is.null(built))
})
