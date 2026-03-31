# ============================================================================== #
# Lipidomics Analysis Pipeline for HEPAFIT
# ============================================================================== #

library(readxl)
library(dplyr)
library(ggplot2)
library(limma)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(mixOmics)
library(ggrepel)
library(reshape2)
library(tidyr)
library(circlize)

setwd("~HEPAFIT/")
dir.create("Results", showWarnings = FALSE)

# ============================================================================== #
# 1. Data Preparation & Cleaning
# ============================================================================== #

# KEGG Annotation
KEGG_annot <- read_excel("Lipidomics.xlsx", sheet = 5) %>%
    t() %>%
    as.data.frame()
colnames(KEGG_annot) <- KEGG_annot[1, ]
KEGG_annot <- KEGG_annot[-1, c(1, 2, 3, 7, 11, 12)]
colnames(KEGG_annot)[4] <- "Individual_Notation"
KEGG_annot$SubclassB <- KEGG_annot$CLASS
KEGG_annot$SubclassB[KEGG_annot$`SUBCLASS A` == "Triacylglycerols"] <- "Triacylglycerols"

# Lipidomics Data
HEPAFIT <- read_excel("Lipidomics.xlsx", sheet = 4) %>% as.data.frame()
HEPAFIT$PATIENT_ID <- sub("^\\d+_", "", HEPAFIT$COD.EXT)

# Filter for subjects with both Baseline and 6 months
ID_availability <- as.data.frame.matrix(table(HEPAFIT$PATIENT_ID, HEPAFIT$Visit))
IDs <- rownames(ID_availability)[ID_availability$Baseline == 1 & ID_availability$`6 months` == 1]
HEPAFIT <- HEPAFIT[HEPAFIT$PATIENT_ID %in% IDs, ]

# Separate and order Timepoints
HEPAFIT_pre <- HEPAFIT[HEPAFIT$Visit == "Baseline", ]
HEPAFIT_pre <- HEPAFIT_pre[order(HEPAFIT_pre$PATIENT_ID), ]

HEPAFIT_post <- HEPAFIT[HEPAFIT$Visit == "6 months", ]
HEPAFIT_post <- HEPAFIT_post[order(HEPAFIT_post$PATIENT_ID), ]

# Metadata
metadata <- read_excel("Lipidomics.xlsx", sheet = 4) %>% as.data.frame()
metadata$PATIENT_ID <- sub("^\\d+_", "", metadata$COD.EXT)
metadata <- metadata[metadata$PATIENT_ID %in% IDs, ]
metadata <- metadata[order(metadata$PATIENT_ID), ]

clinical_meta <- read_excel("Lipidomics.xlsx", sheet = 6) %>% as.data.frame()
clinical_meta$PATIENT_ID <- sub("^\\d+_", "", clinical_meta$COD.EXT)
clinical_meta <- clinical_meta[clinical_meta$PATIENT_ID %in% IDs, ]
clinical_meta <- clinical_meta[order(clinical_meta$PATIENT_ID), ]

Arm <- HEPAFIT_pre$Arm
Sex <- factor(clinical_meta$`SEX (0 Girls, 1 Boys)`, levels = c(0, 1), labels = c("Girls", "Boys"))


# ============================================================================== #
# 2. Normalization Vectors
# ============================================================================== #

# Numeric Matrices
lipid_cols <- 10:270
HEPAFIT_num_pre <- HEPAFIT_pre[, lipid_cols]
HEPAFIT_num_post <- HEPAFIT_post[, lipid_cols]

# Log2 Transform (For true log2FC interpretation in t-tests)
HEPAFIT_log_pre <- log2(HEPAFIT_num_pre)
HEPAFIT_log_post <- log2(HEPAFIT_num_post)
HEPAFIT_log_all <- rbind(HEPAFIT_log_pre, HEPAFIT_log_post)

# Scaled (For PCA / PLS-DA only, doing fold change on scaled data is risky)
HEPAFIT_scaled_all <- scale(HEPAFIT_log_all, center = TRUE, scale = TRUE)
HEPAFIT_scaled_pre <- HEPAFIT_scaled_all[1:length(IDs), ]
HEPAFIT_scaled_post <- HEPAFIT_scaled_all[(length(IDs) + 1):(2 * length(IDs)), ]

# ============================================================================== #
# 3. PCA
# ============================================================================== #

plot_pca <- function(data, group, visit = NULL, title = "PCA") {
    pca_res <- prcomp(data)
    var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)

    mds <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], Group = group)
    if (!is.null(visit)) mds$Visit <- visit

    p <- ggplot(mds, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 1.5, show.legend = TRUE)
    if (!is.null(visit)) p <- p + aes(shape = Visit)

    p + labs(
        title = title,
        x = sprintf("PC1 (%.2f%%)", var_exp[1] * 100),
        y = sprintf("PC2 (%.2f%%)", var_exp[2] * 100)
    ) +
        theme_minimal() +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.2)
}


pdf("Results/PCA Baseline (Scaled)_sex.pdf")
plot_pca(HEPAFIT_scaled_pre, Arm, title = "PCA Baseline by intervention group (Scaled)")
dev.off()

pdf("Results/PCA Baseline sex (Scaled)_sex.pdf")
plot_pca(HEPAFIT_scaled_pre, as.factor(clinical_meta$`SEX (0 Girls, 1 Boys)`), title = "PCA Baseline by sex (Scaled)")
dev.off()


pdf("PCA Post-Exercise (Scaled)_sex.pdf", width = 12, height = 8)
plot_pca(HEPAFIT_scaled_post, Arm, title = "PCA Post-Exercise (Scaled)")
dev.off()

# ============================================================================== #
# 4. PRE vs POST Analysis (Paired Limma per Group)
# ============================================================================== #

run_paired_limma <- function(group_name) {
    # Subset
    idx <- which(Arm == group_name)
    pre <- HEPAFIT_log_pre[idx, ]
    post <- HEPAFIT_log_post[idx, ]

    # Limma requires variables (lipids) as rows, samples as columns
    exprs <- t(rbind(pre, post))

    # Define paired design mapping subjects and timepoints
    subjects <- factor(rep(IDs[idx], 2))
    time <- factor(rep(c("Pre", "Post"), each = length(idx)), levels = c("Pre", "Post"))

    # Paired design matrix: intercept for each subject + time effect
    # Note: Sex is perfectly correlated with `subjects` so it cannot be added as a main effect here.
    # The paired statistical design already naturally corrects for biological sex differences
    # because patients act as their own control.
    design <- model.matrix(~ subjects + time)

    # Fit Limma model
    fit <- lmFit(exprs, design)
    fit <- eBayes(fit, robust = TRUE)

    # Extract results for the 'timePost' coefficient (Post - Pre)
    res <- topTable(fit, coef = "timePost", number = Inf, sort.by = "P", confint = TRUE)
    res$Lipid <- rownames(res)
    res$log2FC <- res$logFC
    res$p_value <- res$P.Value
    res$adj_p_value <- res$adj.P.Val
    res$CI.L <- res$CI.L
    res$CI.R <- res$CI.R
    res$Group <- group_name

    # Merge annot
    res <- merge(res, KEGG_annot, by.x = "Lipid", by.y = "VARIABLE OWL ID")
    res$cat_volcano <- res$SubclassB
    res$cat_volcano[res$p_value > 0.05] <- "ns"
    res$cat_volcano <- factor(res$cat_volcano, levels = c(
        "Glycerolipids", "Glycerophospholipids",
        "Sphingolipids", "Sterols",
        "Triacylglycerols", "ns"
    ))
    return(res)
}

res_paired_ctrl <- run_paired_limma("Control")
res_paired_hipe <- run_paired_limma("HIPE")
res_paired_lipe <- run_paired_limma("LIPE")
res_paired_plus <- run_paired_limma("PLUS")

plot_volcano <- function(df, title) {
    ggplot(df, aes(x = log2FC, y = -log10(p_value), colour = cat_volcano)) +
        geom_vline(xintercept = c(-0.32, 0.32), linetype = "dashed", color = "grey40", size = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.5) +
        geom_point(alpha = 0.85, size = 2) +
        scale_color_manual(
            values = c(
                "Glycerolipids" = "#1f78b4", # Rich Blue
                "Glycerophospholipids" = "#e31a1c", # Deep Red
                "Sphingolipids" = "#33a02c", # Forest Green
                "Sterols" = "#ff7f00", # Vivid Orange
                "Triacylglycerols" = "#6a3d9a", # Royal Purple
                "ns" = "#dce0e5" # Soft background grey
            ), drop = FALSE, name = "Lipid Class"
        ) +
        geom_text_repel(
            data = subset(df, p_value < 0.05 & abs(log2FC) > 0.32),
            aes(label = Individual_Notation),
            size = 3.5,
            box.padding = 0.5,
            point.padding = 0.3,
            max.overlaps = 25,
            segment.color = "grey40",
            segment.size = 0.4,
            min.segment.length = 0,
            show.legend = FALSE,
            fontface = "italic"
        ) +
        labs(
            title = title,
            x = expression(bold("log"[2] * " Fold Change (Post - Pre)")),
            y = expression(bold("-log"[10] * " P-value"))
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.text = element_text(color = "black", size = 9),
            axis.title = element_text(size = 10),
            legend.title = element_text(face = "bold", size = 13),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", size = 1)
        )
}

v_ctrl <- plot_volcano(res_paired_ctrl, "Control: Pre vs Post (Limma)")
v_hipe <- plot_volcano(res_paired_hipe, "HIPE: Pre vs Post (Limma)")
v_lipe <- plot_volcano(res_paired_lipe, "LIPE: Pre vs Post (Limma)")
v_plus <- plot_volcano(res_paired_plus, "PLUS: Pre vs Post (Limma)")

pdf("Results/Volcanos_Pre_Post_Limma_by_Group_sex.pdf", width = 12, height = 8)
arranged_volcanos_prepost <- ggarrange(v_ctrl, v_hipe, v_lipe, v_plus, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
print(arranged_volcanos_prepost)
dev.off()

# Divergent Bar Plots for Significant Lipids (Pre vs Post)
plot_divergent_bars <- function(df, title) {
    # Filter for significant lipids only
    sig_df <- subset(df, p_value < 0.05)

    if (nrow(sig_df) == 0) {
        return(ggplot() +
            theme_void() +
            ggtitle(paste(title, "(No Significant Lipids)")))
    }

    # Calculate counts of up and down regulated lipids per subclass
    counts <- sig_df %>%
        mutate(
            Direction = ifelse(log2FC > 0, "Up", "Down"),
            Count = ifelse(log2FC > 0, 1, -1)
        ) %>%
        group_by(SubclassB, Direction) %>%
        summarise(Total = sum(Count), .groups = "drop") %>%
        mutate(FillGroup = paste(SubclassB, Direction, sep = "_"))

    ggplot(counts, aes(x = Total, y = SubclassB, fill = FillGroup)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(
            values = c(
                "Glycerolipids_Up" = "#1f78b4", # Rich Blue
                "Glycerolipids_Down" = "#a6cee3", # Light Blue
                "Glycerophospholipids_Up" = "#e31a1c", # Deep Red
                "Glycerophospholipids_Down" = "#fb9a99", # Light Red/Pink
                "Sphingolipids_Up" = "#33a02c", # Forest Green
                "Sphingolipids_Down" = "#b2df8a", # Light Green
                "Sterols_Up" = "#ff7f00", # Vivid Orange
                "Sterols_Down" = "#fdbf6f", # Light Orange
                "Triacylglycerols_Up" = "#6a3d9a", # Royal Purple
                "Triacylglycerols_Down" = "#cab2d6" # Light Purple
            )
        ) +
        scale_x_continuous(labels = function(x) sprintf("%.0f", x)) +
        geom_vline(xintercept = 0, color = "black", size = 0.5) +
        theme_minimal(base_size = 10) +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.text.y = element_text(color = "black", size = 9),
            axis.text.x = element_text(color = "black", size = 10),
            plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
            panel.grid.major.y = element_blank()
        ) +
        labs(title = title, x = "Significant Lipids Count (Down / Up)")
}

bar_ctrl <- plot_divergent_bars(res_paired_ctrl, "Control: Pre vs Post")
bar_hipe <- plot_divergent_bars(res_paired_hipe, "HIPE: Pre vs Post")
bar_lipe <- plot_divergent_bars(res_paired_lipe, "LIPE: Pre vs Post")
bar_plus <- plot_divergent_bars(res_paired_plus, "PLUS: Pre vs Post")

all_paired_results <- rbind(res_paired_ctrl, res_paired_hipe, res_paired_lipe, res_paired_plus)
bar_summary <- plot_divergent_bars(all_paired_results, "Significant Lipids per Subclass (p < 0.05)")

pdf("Results/Divergent_Bars_Pre_Post_Limma_Summary_sex.pdf", width = 8, height = 4)
print(bar_summary)
dev.off()

write.csv(res_paired_ctrl, "Results/DEA_paired_ctrl_sex.csv", row.names = FALSE)
write.csv(res_paired_hipe, "Results/DEA_paired_hipe_sex.csv", row.names = FALSE)
write.csv(res_paired_lipe, "Results/DEA_paired_lipe_sex.csv", row.names = FALSE)
write.csv(res_paired_plus, "Results/DEA_paired_plus_sex.csv", row.names = FALSE)

# ============================================================================== #
# 5. POST-EXERCISE Analysis (Interventions vs Control using Limma)
# ============================================================================== #

run_post_limma <- function(group_name) {
    idx_group <- which(Arm == group_name)
    idx_ctrl <- which(Arm == "Control")

    grp <- HEPAFIT_log_post[idx_group, ]
    ctrl <- HEPAFIT_log_post[idx_ctrl, ]

    # Limma requires variables as rows
    exprs <- t(rbind(ctrl, grp))

    # Define group factor (Control as reference)
    condition <- factor(c(rep("Control", length(idx_ctrl)), rep(group_name, length(idx_group))),
        levels = c("Control", group_name)
    )

    # Bring in Sex for the matched patients to use as a model covariate
    sex_group <- Sex[c(idx_ctrl, idx_group)]

    # Unpaired design matrix controlling for Sex (if there is variance in both groups to avoid model singularity)
    if (length(levels(sex_group)[table(sex_group) > 0]) > 1) {
        design <- model.matrix(~ condition + sex_group)
    } else {
        design <- model.matrix(~condition)
    }

    fit <- lmFit(exprs, design)
    fit <- eBayes(fit, robust = TRUE)

    # Coef is condition[group_name]
    coef_name <- paste0("condition", group_name)
    res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", confint = TRUE)
    res$Lipid <- rownames(res)
    res$log2FC <- res$logFC
    res$p_value <- res$P.Value
    res$adj_p_value <- res$adj.P.Val
    res$CI.L <- res$CI.L
    res$CI.R <- res$CI.R
    res$Group <- group_name

    # Merge annot
    res <- merge(res, KEGG_annot, by.x = "Lipid", by.y = "VARIABLE OWL ID")
    res$cat_volcano <- res$SubclassB
    res$cat_volcano[res$p_value > 0.05] <- "ns"
    res$cat_volcano <- factor(res$cat_volcano, levels = c(
        "Glycerolipids", "Glycerophospholipids",
        "Sphingolipids", "Sterols",
        "Triacylglycerols", "ns"
    ))
    return(res)
}

res_post_hipe <- run_post_limma("HIPE")
res_post_lipe <- run_post_limma("LIPE")
res_post_plus <- run_post_limma("PLUS")


write.csv(res_post_hipe, "Results/DEA_post_hipe_sex.csv", row.names = FALSE)
write.csv(res_post_lipe, "Results/DEA_post_lipe_sex.csv", row.names = FALSE)
write.csv(res_post_plus, "Results/DEA_post_plus_sex.csv", row.names = FALSE)


all_post_results <- rbind(res_post_hipe, res_post_lipe, res_post_plus)

# Volcanos for Post vs Control
vp_hipe <- plot_volcano(res_post_hipe, "HIPE vs Control (Post)")
vp_lipe <- plot_volcano(res_post_lipe, "LIPE vs Control (Post)")
vp_plus <- plot_volcano(res_post_plus, "PLUS vs Control (Post)")


# Significant log2FC
sig_lipids <- unique(all_post_results$Individual_Notation[all_post_results$p_value < 0.05])
df_sig_post <- all_post_results %>% filter(Individual_Notation %in% sig_lipids)

# Circular Plot for Significant Results
# Calculate significance stars to distinguish between the groups!
df_sig_post$sig_star <- ifelse(df_sig_post$p_value < 0.001, "***",
    ifelse(df_sig_post$p_value < 0.01, "**",
        ifelse(df_sig_post$p_value < 0.05, "*", "")
    )
)

# Order Individual_Notation by Class
lipid_order <- df_sig_post %>%
    group_by(Individual_Notation) %>%
    summarise(SubclassB = first(SubclassB), sort_val = max(abs(log2FC))) %>%
    arrange(SubclassB, -sort_val)

df_sig_post$Individual_Notation <- factor(df_sig_post$Individual_Notation, levels = lipid_order$Individual_Notation)

# Determine the Y position for the outer label ring (added extra padding for the asterisks)
ring_y <- max(df_sig_post$CI.R, na.rm = TRUE) + 0.6

# Calculate angles for labels to radiate gracefully outward
lipid_order$id <- seq(1, nrow(lipid_order))
lipid_order$angle <- 90 - 360 * lipid_order$id / nrow(lipid_order)
lipid_order$hjust <- ifelse(lipid_order$angle < -90, 1, 0)
lipid_order$angle <- ifelse(lipid_order$angle < -90, lipid_order$angle + 180, lipid_order$angle)

circular_plot <- ggplot(df_sig_post, aes(x = Individual_Notation, y = log2FC, color = SubclassB)) +
    # Add baseline 0 reference
    geom_hline(yintercept = 0, color = "grey30", linetype = "dashed", size = 0.8) +
    # Expand Y axis to fit our labels using a scale modifier
    scale_y_continuous(limits = c(min(df_sig_post$CI.L, na.rm = TRUE) - 0.4, ring_y + 0.7)) +
    # Plot confidence intervals as lines and dots as estimates (dodged so groups don't overlap)
    # We fade non-significant points so the reader only focuses on the significant ones
    geom_linerange(aes(ymin = CI.L, ymax = CI.R, group = Group), position = position_dodge(width = 0.7), alpha = ifelse(df_sig_post$p_value < 0.05, 0.6, 0.2), size = 0.6) +
    geom_point(aes(shape = Group, group = Group), position = position_dodge(width = 0.7), size = 3, alpha = ifelse(df_sig_post$p_value < 0.05, 0.9, 0.3)) +
    # Add significance asterisks at the edge of the confidence intervals
    geom_text(aes(label = sig_star, y = ifelse(log2FC > 0, CI.R + 0.15, CI.L - 0.15), group = Group), position = position_dodge(width = 0.7), size = 4, color = "black", vjust = 0.7) +
    # Add radiating explicit label text for each lipid on the outer ring
    geom_text(
        data = lipid_order,
        aes(x = Individual_Notation, y = ring_y, label = Individual_Notation, angle = angle, hjust = hjust),
        size = 3.2, color = "black", fontface = "italic", inherit.aes = FALSE
    ) +
    # Consistent color scale with previously updated Volcano standard
    scale_color_manual(
        values = c(
            "Glycerolipids" = "#1f78b4", # Rich Blue
            "Glycerophospholipids" = "#e31a1c", # Deep Red
            "Sphingolipids" = "#33a02c", # Forest Green
            "Sterols" = "#ff7f00", # Vivid Orange
            "Triacylglycerols" = "#6a3d9a" # Royal Purple
        ), name = "Lipid Class"
    ) +
    # Polar coordinate mapping
    coord_polar() +
    theme_minimal(base_size = 10) +
    theme(
        axis.text.x = element_blank(), # Hiding exact x-axis text to prevent horrible overlapping natively
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major.x = element_line(color = "grey80", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
        legend.position = "right",
        legend.box = "vertical"
    ) +
    labs(
        title = "Post-Exercise vs Control: Circular log2FC Estimates",
        x = "",
        y = expression(bold("log"[2] * " Fold Change (95% CI)"))
    )

pdf("Results/Post_Exercise_vs_Control_Limma_Circular_Plot_sex.pdf", width = 11, height = 9)
print(circular_plot)
dev.off()


# ============================================================================== #
# 6. Sparse PLS-DA (sPLS-DA)
# ============================================================================== #

# Function to run sPLS-DA and immediately return the individual score plot
run_splsda_indiv <- function(X, Y, title_name, colors = NULL, multilevel = NULL) {
    # Run sPLS-DA (keeps the top 20 features per component)
    if (!is.null(multilevel)) {
        model <- splsda(X, Y, ncomp = 2, keepX = c(20, 20), multilevel = multilevel)
    } else {
        model <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
    }

    # Construct the plot string dynamically to accommodate custom colors natively
    if (!is.null(colors)) {
        p_indiv <- plotIndiv(model,
            comp = c(1, 2),
            group = Y,
            legend = TRUE,
            ellipse = TRUE,
            title = paste("PLS-DA", title_name),
            pch = 16,
            col.per.group = colors,
            size.title = 14
        )
    } else {
        p_indiv <- plotIndiv(model,
            comp = c(1, 2),
            group = Y,
            legend = TRUE,
            ellipse = TRUE,
            title = paste("PLS-DA", title_name),
            pch = 16,
            size.title = 14
        )
    }

    # Overwrite the native mixOmics layout with a beautiful publication-ready ggplot theme
    custom_plot <- p_indiv$graph +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            axis.text = element_text(color = "black", size = 9),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", size = 1)
        )

    # Return the customized ggplot object natively to render in R directly
    return(custom_plot)
}

# --- PRE vs POST (4 separate plots for each group)
run_splsda_per_arm <- function(group_name) {
    idx <- which(Arm == group_name)
    # Subset scaled data for 'pre' and 'post' for only this group
    pre <- HEPAFIT_scaled_all[idx, ]
    post <- HEPAFIT_scaled_all[idx + length(Arm), ]
    X_grp <- rbind(pre, post)

    Y_grp <- factor(rep(c("Baseline", "Post-Exercise"), each = length(idx)), levels = c("Baseline", "Post-Exercise"))
    patient_grp <- data.frame(sample = rep(IDs[idx], 2))

    # Run sPLS-DA using paired approach and matched colors (Blue and Orange)
    plot_arm <- run_splsda_indiv(X_grp, Y_grp, title_name = group_name, colors = c("#2c7bb6", "#ff7f00"), multilevel = patient_grp)
    return(plot_arm)
}

# Generate plots for each group
plsda_ctrl <- run_splsda_per_arm("Control")
plsda_hipe <- run_splsda_per_arm("HIPE")
plsda_lipe <- run_splsda_per_arm("LIPE")
plsda_plus <- run_splsda_per_arm("PLUS")

# Visualize directly in R (2x2 grid)
arranged_prepost_splsda <- ggarrange(plsda_ctrl, plsda_hipe, plsda_lipe, plsda_plus, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
print(arranged_prepost_splsda)

pdf("Results/arranged_prepost_splsda_sex.pdf")
print(arranged_prepost_splsda)
dev.off()



# ============================================================================== #
# 7. Intervention Target Overlap & Top Lipid Boxplots
# ============================================================================== #

# Top 10 Lipids Boxplots - Strategy 2 (Post-Exercise vs Control)
# Find the exact top 10 most significant lipids comparing interventions vs control
top_post_lipids <- all_post_results %>%
    arrange(p_value) %>%
    pull(Individual_Notation) %>%
    unique() %>%
    head(10)

boxplot_data_post <- data.frame()
for (lipid_name in top_post_lipids) {
    lipid_row <- KEGG_annot[KEGG_annot$Individual_Notation == lipid_name, ]
    if (nrow(lipid_row) > 0) {
        lipid_id <- lipid_row$`VARIABLE OWL ID`[1]

        df <- data.frame(
            Patient = IDs,
            Arm = Arm,
            Value = HEPAFIT_log_post[[lipid_id]],
            Lipid = lipid_name
        )
        boxplot_data_post <- rbind(boxplot_data_post, df)
    }
}

boxplot_data_post$Arm <- factor(boxplot_data_post$Arm, levels = c("Control", "HIPE", "LIPE", "PLUS"))

p_box_post <- ggplot(boxplot_data_post, aes(x = Arm, y = Value, fill = Arm)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(color = "black", size = 1, alpha = 0.6, width = 0.15) +
    facet_wrap(~Lipid, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("Control" = "grey", "HIPE" = "#F8766D", "LIPE" = "#00BA38", "PLUS" = "#619CFF")) +
    theme_bw(base_size = 9) +
    labs(title = "Top 10 Discriminative Lipids: Interventions vs Control (Post)", x = "", y = "log2 Intensity") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
    )

pdf("Results/Boxplots_Top10_PostControl_sex.pdf", width = 14, height = 7)
print(p_box_post)
dev.off()

# ============================================================================== #
# 8. Clinical Correlation Heatmaps (Lipid-to-Phenotype Trajectories)
# ============================================================================== #

# Load the clinical data explicitly mapped on sheet 6
clinical_data <- read_excel("Lipidomics.xlsx", sheet = 6) %>% as.data.frame()
clinical_data$PATIENT_ID <- sub("^\\d+_", "", clinical_data$COD.EXT)

# Merge so patient order matches our main 'IDs' vector perfectly
clinical_matched <- data.frame(PATIENT_ID = IDs) %>%
    left_join(clinical_data, by = "PATIENT_ID")

# Calculate Spearman correlations for BOTH sets of Top 10 Lipids combined
top_cor_lipids <- unique(c(top_prepost_lipids, top_post_lipids))

# We specifically want delta values where available, and baseline for others
# Excluded categorical text/binary variables ("NAFLD", "DIETA", "CODE IR") in favor of mathematically robust, continuous variables for Spearman correlation
clinical_vars_to_keep <- c("DIF_homa", "DIF_CAP", "BMI", "BODY_FAT")
delta_clinical <- clinical_matched[, clinical_vars_to_keep]
rownames(delta_clinical) <- IDs
colnames(delta_clinical) <- c("Δ HOMA-IR", "Δ CAP Liver Fat", "Baseline BMI", "Baseline Body Fat")

delta_lipid_matrix <- matrix(NA, nrow = length(IDs), ncol = length(top_cor_lipids))
rownames(delta_lipid_matrix) <- IDs
colnames(delta_lipid_matrix) <- top_cor_lipids

# Loop and build log2FC (Post - Pre) for each patient
for (i in seq_along(top_cor_lipids)) {
    lipid_name <- top_cor_lipids[i]
    lipid_row <- KEGG_annot[KEGG_annot$Individual_Notation == lipid_name, ]
    if (nrow(lipid_row) > 0) {
        lipid_id <- lipid_row$`VARIABLE OWL ID`[1]

        # log2FC is inherently Post - Pre for log transformed data
        delta_lipid_matrix[, i] <- HEPAFIT_log_post[[lipid_id]] - HEPAFIT_log_pre[[lipid_id]]
    }
}

# Clean overlapping missing values
valid_patients <- rowSums(!is.na(delta_clinical)) > 3
delta_clinical_clean <- delta_clinical[valid_patients, ]
delta_lipid_clean <- delta_lipid_matrix[valid_patients, ]

# Calculate Spearman Correlation Matrix
cor_matrix <- cor(delta_lipid_clean, delta_clinical_clean, method = "spearman", use = "pairwise.complete.obs")

# Compute Exact P-Values to overlay significance asterisks
p_matrix <- matrix(NA, nrow = ncol(delta_lipid_clean), ncol = ncol(delta_clinical_clean))
rownames(p_matrix) <- colnames(delta_lipid_clean)
colnames(p_matrix) <- colnames(delta_clinical_clean)

for (i in seq_len(ncol(delta_lipid_clean))) {
    for (j in seq_len(ncol(delta_clinical_clean))) {
        valid_idx <- !is.na(delta_lipid_clean[, i]) & !is.na(delta_clinical_clean[, j])
        if (sum(valid_idx) > 3) {
            test_res <- cor.test(delta_lipid_clean[valid_idx, i], delta_clinical_clean[valid_idx, j], method = "spearman", exact = FALSE)
            p_matrix[i, j] <- test_res$p.value
        } else {
            p_matrix[i, j] <- 1
        }
    }
}

# Create a single concatenated heatmap partitioned by group
groups_to_plot <- c("Control", "HIPE", "LIPE", "PLUS")

cor_matrix_combined <- NULL
p_matrix_combined <- NULL
group_annotation_vector <- c()
colnames_combined <- c()

for (group_name in groups_to_plot) {
    # Isolate patients belonging to this group
    idx_group <- which(Arm[valid_patients] == group_name)

    if (length(idx_group) > 3) { # Need minimum samples for spearman correlation
        group_lipid <- delta_lipid_clean[idx_group, ]
        group_clinical <- delta_clinical_clean[idx_group, ]

        # Recalculate Correlation Matrix
        cor_matrix_grp <- cor(group_lipid, group_clinical, method = "spearman", use = "pairwise.complete.obs")

        # Recalculate P-Values
        p_matrix_grp <- matrix(NA, nrow = ncol(group_lipid), ncol = ncol(group_clinical))
        rownames(p_matrix_grp) <- colnames(group_lipid)
        colnames(p_matrix_grp) <- colnames(group_clinical)

        for (i in seq_len(ncol(group_lipid))) {
            for (j in seq_len(ncol(group_clinical))) {
                valid_idx_grp <- !is.na(group_lipid[, i]) & !is.na(group_clinical[, j])
                if (sum(valid_idx_grp) > 3) {
                    test_res <- cor.test(group_lipid[valid_idx_grp, i], group_clinical[valid_idx_grp, j], method = "spearman", exact = FALSE)
                    p_matrix_grp[i, j] <- test_res$p.value
                } else {
                    p_matrix_grp[i, j] <- 1
                }
            }
        }

        # Standardize missing correlation values
        cor_matrix_grp[is.na(cor_matrix_grp)] <- 0

        # Append to our combined matrices
        if (is.null(cor_matrix_combined)) {
            cor_matrix_combined <- cor_matrix_grp
            p_matrix_combined <- p_matrix_grp
        } else {
            cor_matrix_combined <- cbind(cor_matrix_combined, cor_matrix_grp)
            p_matrix_combined <- cbind(p_matrix_combined, p_matrix_grp)
        }

        group_annotation_vector <- c(group_annotation_vector, rep(group_name, ncol(group_clinical)))
        colnames_combined <- c(colnames_combined, colnames(group_clinical))
    }
}

if (!is.null(cor_matrix_combined)) {
    colnames(cor_matrix_combined) <- colnames_combined

    # Define an annotation bar mapping the columns to their respective groups
    col_annot <- HeatmapAnnotation(
        Group = group_annotation_vector,
        col = list(Group = c("Control" = "grey", "HIPE" = "#F8766D", "LIPE" = "#00BA38", "PLUS" = "#619CFF")),
        annotation_name_side = "left"
    )

    pdf("Results/Clinical_Correlations_Heatmap_All_Groups_sex.pdf", width = 14, height = 8)

    tryCatch(
        {
            ht_clinical_post <- Heatmap(cor_matrix_combined,
                name = "Spearman Rho\nCorrelation",
                col = col_fun_cor,
                top_annotation = col_annot,
                column_split = factor(group_annotation_vector, levels = groups_to_plot),
                column_title = "Clinical Associations segmented by Intervention Group",
                row_title = "Δ log2FC (Most Significant Lipids)",
                column_names_rot = 45,
                column_names_centered = FALSE,
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 9, fontface = "bold"),
                rect_gp = gpar(col = "white", lwd = 1),
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                cell_fun = function(j, i, x, y, width, height, fill) {
                    if (!is.na(p_matrix_combined[i, j]) && p_matrix_combined[i, j] < 0.05) {
                        grid.text("*", x, y, vjust = 0.7, gp = gpar(fontsize = 9, col = "black", fontface = "bold"))
                    }
                }
            )
            print(ht_clinical_post)
        },
        error = function(e) {
            cat("Could not generate merged group heatmap - likely missing variance.\n")
        }
    )

    dev.off()
}

# ============================================================================== #
# 9. Lipid Class Enrichment Analysis (Hypergeometric Test)
# ============================================================================== #

# Calculate Background: Total number of lipids in each subclass measured in the entire dataset
background_subclasses <- KEGG_annot %>%
    filter(!is.na(SubclassB)) %>%
    group_by(SubclassB) %>%
    summarise(Total_Measured = n()) %>%
    filter(SubclassB %in% c("Glycerolipids", "Glycerophospholipids", "Sphingolipids", "Sterols", "Triacylglycerols"))

# Total lipids matched to SubclassB
N_total <- sum(background_subclasses$Total_Measured)

# We define a function to compute hypergeometric p-values for a given set of significant results
calculate_enrichment <- function(results_df, group_name) {
    # 1. Total Significant Lipids for this group (k)
    k_total_sig <- sum(results_df$p_value < 0.05, na.rm = TRUE)

    if (k_total_sig == 0) {
        return(NULL)
    } # No significant lipids, no enrichment possible

    # 2. Count significant lipids per subclass (x)
    sig_subclasses <- results_df %>%
        filter(p_value < 0.05) %>%
        group_by(SubclassB) %>%
        summarise(Sig_Count = n())

    # Merge with background dataset
    enrichment <- left_join(background_subclasses, sig_subclasses, by = "SubclassB")
    enrichment$Sig_Count[is.na(enrichment$Sig_Count)] <- 0

    # Calculate hypergeometric p-value for each class
    # x: number of significant lipids in class
    # m: total measured lipids in class
    # n: total measured lipids NOT in class (N_total - m)
    # k: total significant lipids overall

    enrichment$p_hyper <- apply(enrichment, 1, function(row) {
        x <- as.numeric(row["Sig_Count"])
        m <- as.numeric(row["Total_Measured"])
        n <- N_total - m
        k <- k_total_sig

        # phyper calculates P(X <= x). We want P(X >= x) for enrichment.
        # So we use phyper(x - 1, m, n, k, lower.tail = FALSE)
        p_val <- phyper(x - 1, m, n, k, lower.tail = FALSE)
        return(p_val)
    })

    # FDR correction
    enrichment$FDR <- p.adjust(enrichment$p_hyper, method = "BH")

    # Calculate Fold Enrichment = (x/k) / (m/N)
    enrichment$Fold_Enrichment <- (enrichment$Sig_Count / k_total_sig) / (enrichment$Total_Measured / N_total)

    enrichment$Group <- group_name
    return(enrichment)
}

# Run the enrichment across the Pre vs Post paired analyses for each group
enrich_ctrl <- calculate_enrichment(res_paired_ctrl, "Control")
enrich_hipe <- calculate_enrichment(res_paired_hipe, "HIPE")
enrich_lipe <- calculate_enrichment(res_paired_lipe, "LIPE")
enrich_plus <- calculate_enrichment(res_paired_plus, "PLUS")

# Combine all results
all_enrichment <- bind_rows(enrich_ctrl, enrich_hipe, enrich_lipe, enrich_plus)

# Filter out classes with no significant hits to avoid plotting empty bubbles
all_enrichment <- all_enrichment %>% filter(Sig_Count > 0)
all_enrichment$Group <- factor(all_enrichment$Group, levels = c("Control", "HIPE", "LIPE", "PLUS"))

# Build a Bubble Plot for Enrichment (Pre vs Post)
p_enrichment <- ggplot(all_enrichment, aes(x = Group, y = SubclassB)) +
    geom_point(aes(size = Fold_Enrichment, fill = -log10(p_hyper)), shape = 21, color = "black", stroke = 1) +
    scale_fill_gradient(low = "white", high = "#e31a1c", name = "-log10(p-value)") +
    scale_size_continuous(range = c(5, 20), limits = c(0, NA), name = "Fold Enrichment ratio") +
    theme_bw(base_size = 10) +
    labs(
        title = "Structural Over-Representation Analysis",
        subtitle = "Hypergeometric Lipid Class Enrichment (Pre vs Post)",
        x = "Group",
        y = "Lipid Class"
    ) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic", margin = margin(b = 15)),
        axis.text.x = element_text(face = "bold", size = 9, color = "grey30"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title = element_text(size = 10),
        panel.grid.major.x = element_line(linetype = "dotted", color = "grey60", size = 0.8),
        panel.grid.major.y = element_line(linetype = "dotted", color = "grey60", size = 0.8),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 10)
    )

pdf("Results/Lipid_Class_Enrichment_Pre_vs_Post_BubblePlot_sex.pdf", width = 8, height = 6.5)
print(p_enrichment)
dev.off()


# Run the enrichment across the Post vs Control analyses for each intervention
enrich_post_hipe <- calculate_enrichment(res_post_hipe, "HIPE")
enrich_post_lipe <- calculate_enrichment(res_post_lipe, "LIPE")
enrich_post_plus <- calculate_enrichment(res_post_plus, "PLUS")

all_enrichment_post <- bind_rows(enrich_post_hipe, enrich_post_lipe, enrich_post_plus)

if (nrow(all_enrichment_post) > 0) {
    all_enrichment_post <- all_enrichment_post %>% filter(Sig_Count > 0)

    if (nrow(all_enrichment_post) > 0) {
        all_enrichment_post$Group <- factor(all_enrichment_post$Group, levels = c("HIPE", "LIPE", "PLUS"))

        p_enrichment_post <- ggplot(all_enrichment_post, aes(x = Group, y = SubclassB)) +
            geom_point(aes(size = Fold_Enrichment, fill = -log10(p_hyper)), shape = 21, color = "black", stroke = 1) +
            scale_fill_gradient(low = "white", high = "#e31a1c", name = "-log10(p-value)") +
            scale_size_continuous(range = c(5, 20), limits = c(0, NA), name = "Fold Enrichment ratio") +
            theme_bw(base_size = 10) +
            labs(
                title = "Structural Over-Representation Analysis",
                subtitle = "Hypergeometric Lipid Class Enrichment (Post vs Control)",
                x = "Exercise Intervention vs Control",
                y = "Lipid Class"
            ) +
            theme(
                plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
                plot.subtitle = element_text(hjust = 0.5, size = 9, face = "italic", margin = margin(b = 15)),
                axis.text.x = element_text(face = "bold", size = 9, color = "grey30"),
                axis.text.y = element_text(size = 9, face = "bold", color = "black"),
                axis.title = element_text(size = 10),
                panel.grid.major.x = element_line(linetype = "dotted", color = "grey60", size = 0.8),
                panel.grid.major.y = element_line(linetype = "dotted", color = "grey60", size = 0.8),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(color = "black", size = 1, fill = NA),
                legend.title = element_text(size = 9, face = "bold"),
                legend.text = element_text(size = 10),
                legend.position = "bottom",
                legend.box = "vertical"
            )

        pdf("Results/Lipid_Class_Enrichment_Post_vs_Control_BubblePlot_sex.pdf", width = 8, height = 7)
        print(p_enrichment_post)
        dev.off()
    }
}


# ============================================================================== #
# 10. Export for LION/web Ontology Enrichment
# ============================================================================== #
# LION/web requires a "Target List" (significant lipids) and a "Background/Universe"
# We will create simple CSV files with the pure lipid names that match LION nomenclature.

# Background universe: All measured lipids that have an individual notation mapping
background_lipids <- unique(KEGG_annot$Individual_Notation[!is.na(KEGG_annot$Individual_Notation)])
write.table(background_lipids, "Results/LION_Background_Universe_sex.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Target list: HIPE
write.table(sig_hipe, "Results/LION_Target_HIPE_sex.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Target list: LIPE
write.table(sig_lipe, "Results/LION_Target_LIPE_sex.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Target list: PLUS
write.table(sig_plus, "Results/LION_Target_PLUS_sex.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Target list: Post-Exercise vs Control combined
sig_post <- all_post_results %>%
    filter(p_value < 0.05) %>%
    pull(Individual_Notation) %>%
    unique()
write.table(sig_post, "Results/LION_Target_Post_vs_Control_sex.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("\nLION/web export complete. Upload these CSVs directly to www.lipidontology.com\n")

# ============================================================================== #
# 11. Save Differential Expression Analysis (DEA) Results
# ============================================================================== #

write.csv(all_paired_results, "Results/DEA_Pre_vs_Post_sex.csv", row.names = FALSE)
write.csv(all_post_results, "Results/DEA_Post_vs_Control_sex.csv", row.names = FALSE)

cat("DEA results successfully saved to the Results/ directory.\n")

# ============================================================================== #
# 12. Baseline Sex Differences (Supplementary)
# ============================================================================== #
# Checking if there are global baseline differences between Boys and Girls.
cat("\nRunning isolated supplementary baseline analysis for Boys vs Girls...\n")

# Use the log-transformed baseline data matrices
exprs_base <- t(HEPAFIT_log_pre)

# Verify we have variation in Sex to proceed
if (length(levels(Sex)[table(Sex) > 0]) > 1) {
    # Design matrix for Sex at baseline (Reference: Girls)
    design_sex <- model.matrix(~Sex)

    fit_sex <- lmFit(exprs_base, design_sex)
    fit_sex <- eBayes(fit_sex, robust = TRUE)

    # Extract results for Boys vs Girls (coefficient implicitly named SexBoys as Girls is the first factor level)
    res_sex <- topTable(fit_sex, coef = "SexBoys", number = Inf, sort.by = "P", confint = TRUE)
    res_sex$Lipid <- rownames(res_sex)
    res_sex$log2FC <- res_sex$logFC
    res_sex$p_value <- res_sex$P.Value
    res_sex$adj_p_value <- res_sex$adj.P.Val
    res_sex$CI.L <- res_sex$CI.L
    res_sex$CI.R <- res_sex$CI.R

    # Merge KEGG annot
    res_sex <- merge(res_sex, KEGG_annot, by.x = "Lipid", by.y = "VARIABLE OWL ID")
    res_sex$cat_volcano <- res_sex$SubclassB
    res_sex$cat_volcano[res_sex$p_value > 0.05] <- "ns"
    res_sex$cat_volcano <- factor(res_sex$cat_volcano, levels = c(
        "Glycerolipids", "Glycerophospholipids",
        "Sphingolipids", "Sterols",
        "Triacylglycerols", "ns"
    ))

    # Save the table
    write.csv(res_sex, "Results/DEA_Baseline_Sex_Differences_sex.csv", row.names = FALSE)

    # Generate a volcano plot
    vp_sex <- plot_volcano(res_sex, "Baseline Differences: Boys vs Girls")

    pdf("Results/Volcano_Baseline_Sex_Differences_sex.pdf", width = 7, height = 5)
    print(vp_sex)
    dev.off()

    # --- Clinical Baseline Characteristics by Sex ---
    cat("Calculating clinical baseline characteristics by Sex...\n")

    # Target continuous clinical variables measured at baseline
    clinical_vars_cont <- c("Age", "BMI", "BODY_FAT", "HOMA_IR", "CAP_PRE")
    # Target categorical clinical variables measured at baseline
    clinical_vars_cat <- c("HOMA_cat>2.6", "CODE IR", "NAFLD Higado graso (CAP>249)", "NAFLD 2 (CAP>225)")

    clinical_res <- data.frame(
        Variable = character(),
        Girls_Stat = character(),
        Boys_Stat = character(),
        P_Value = numeric(),
        stringsAsFactors = FALSE
    )

    # Calculate group descriptive statistics and perform an independent t-test for continuous vars
    for (var in clinical_vars_cont) {
        if (var %in% colnames(clinical_meta)) {
            val_girls <- as.numeric(clinical_meta[[var]][Sex == "Girls"])
            val_boys <- as.numeric(clinical_meta[[var]][Sex == "Boys"])

            stat_girls <- sprintf("%.2f (%.2f)", mean(val_girls, na.rm = TRUE), sd(val_girls, na.rm = TRUE))
            stat_boys <- sprintf("%.2f (%.2f)", mean(val_boys, na.rm = TRUE), sd(val_boys, na.rm = TRUE))

            # Validate enough complete observations to perform a t-test securely
            if (sum(!is.na(val_girls)) > 2 && sum(!is.na(val_boys)) > 2) {
                t_test <- t.test(val_girls, val_boys)
                p_val <- t_test$p.value
            } else {
                p_val <- NA
            }

            clinical_res <- rbind(clinical_res, list(
                Variable = var,
                Girls_Stat = paste0(stat_girls, " [Mean (SD)]"),
                Boys_Stat = paste0(stat_boys, " [Mean (SD)]"),
                P_Value = p_val
            ))
        }
    }

    # Calculate counts/percentages and perform Fisher's exact test for categorical vars
    for (var in clinical_vars_cat) {
        if (var %in% colnames(clinical_meta)) {
            val_all <- clinical_meta[[var]]
            valid_idx <- !is.na(val_all)

            if (sum(valid_idx) > 0) {
                tbl <- table(val_all[valid_idx], Sex[valid_idx])

                if (nrow(tbl) > 1 && ncol(tbl) > 1) {
                    # Suppress warnings from fisher in small samples
                    suppressWarnings({
                        f_test <- fisher.test(tbl)
                        p_val <- f_test$p.value
                    })
                } else {
                    p_val <- NA
                }

                # Find the positive class (usually "1" or "YES" or "TRUE"). We pick the last factor level
                lvl_to_report <- rownames(tbl)[length(rownames(tbl))]

                count_girls <- sum(val_all[Sex == "Girls"] == lvl_to_report, na.rm = TRUE)
                total_girls <- sum(!is.na(val_all[Sex == "Girls"]))
                perc_girls <- ifelse(total_girls > 0, count_girls / total_girls * 100, 0)

                count_boys <- sum(val_all[Sex == "Boys"] == lvl_to_report, na.rm = TRUE)
                total_boys <- sum(!is.na(val_all[Sex == "Boys"]))
                perc_boys <- ifelse(total_boys > 0, count_boys / total_boys * 100, 0)

                stat_girls <- sprintf("%d (%.1f%%)", count_girls, perc_girls)
                stat_boys <- sprintf("%d (%.1f%%)", count_boys, perc_boys)

                clinical_res <- rbind(clinical_res, list(
                    Variable = paste0(var, " (=", lvl_to_report, ")"),
                    Girls_Stat = paste0(stat_girls, " [n (%)]"),
                    Boys_Stat = paste0(stat_boys, " [n (%)]"),
                    P_Value = p_val
                ))
            }
        }
    }

    # Round P Value for cleaner export
    clinical_res$P_Value <- round(clinical_res$P_Value, 4)
    write.csv(clinical_res, "Results/Clinical_Baseline_Sex_Differences_sex.csv", row.names = FALSE)

    cat("Baseline Sex Differences analysis successfully completed and exported.\n")
} else {
    cat("Skipping baseline differences: Not enough variation in the Sex variable.\n")
}
# ============================================================================== #
# 13. Master Panel Pre vs Post
# ============================================================================== #
cat("\nGenerating absolute master layout panel for all Pre vs Post comparisons...\n")

# Combine the 4 primary visualization assets
# A) arranged_prepost_splsda (2x2 grid)
# B) arranged_volcanos_prepost (2x2 grid)
# C) bar_summary (1x1 wide plot)
# D) p_enrichment (1x1 wide bubble plot)

# To balance sizes, top row gets the multi-panels, bottom row gets the single plots
# We use ggarrange to pack them together with appropriate heights matching the visual density
master_prepost <- ggarrange(
    arranged_prepost_splsda,
    arranged_volcanos_prepost,
    bar_summary,
    p_enrichment,
    labels = c("A", "B", "C", "D"),
    ncol = 2, nrow = 2,
    heights = c(1.5, 1) # Give the 2x2 grids more vertical space
)

pdf("Results/Master_Panel_Pre_vs_Post_sex.pdf", width = 18, height = 15)
print(master_prepost)
dev.off()

cat("Master panel generated successfully.\n")

# ============================================================================== #
# 14. Master Panel Post vs Control
# ============================================================================== #
cat("\nGenerating absolute master layout panel for all Post vs Control comparisons...\n")

# Need to explicitly capture the ComplexHeatmap as a grid graphic object for ggarrange, placing legends at the bottom
ht_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(ht_clinical_post, heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))

# A) arranged_volcanos_post (3x1 grid)
# B) p_enrichment_post (bubble)
# C) circular_plot (polar coord)
# D) ht_grob (Heatmap)

# Top row gets the wide volcano 3-grid and the bubble plot
# Bottom row gets the tall circular plot and the complex heatmap
master_post_top <- ggarrange(
    arranged_volcanos_post,
    p_enrichment_post,
    labels = c("A", "B"),
    ncol = 2, nrow = 1,
    widths = c(1.5, 1) # Volcano grid is natively wider (3 vs 1 columns)
)

master_post_bottom <- ggarrange(
    circular_plot,
    ht_grob,
    labels = c("C", "D"),
    ncol = 2, nrow = 1,
    widths = c(1.2, 1) # Give circular plot more width horizontally to make it bigger
)

master_post <- ggarrange(
    master_post_top,
    master_post_bottom,
    ncol = 1, nrow = 2,
    heights = c(1, 1.5) # Give the bottom row (circular plot + heatmap) more vertical space
)

pdf("Results/Master_Panel_Post_vs_Control_sex.pdf", width = 19, height = 18)
print(master_post)
dev.off()

cat("Master panel exactly for post comparisons generated successfully.\n")
# ============================================================================== #
# 15. Baseline Differences Among Groups
# ============================================================================== #
cat("\nRunning baseline differences analysis among groups (Control, HIPE, LIPE, PLUS)...\n")

# --- 14.1 Clinical Baseline Characteristics Among Groups ---
clinical_res_group <- data.frame(
    Variable = character(),
    Control = character(),
    HIPE = character(),
    LIPE = character(),
    PLUS = character(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
)

# Function to get formatted string for continuous variables
get_cont_stat <- function(var, group) {
    val <- as.numeric(clinical_meta[[var]][Arm == group])
    sprintf("%.2f (%.2f)", mean(val, na.rm = TRUE), sd(val, na.rm = TRUE))
}

for (var in clinical_vars_cont) {
    if (var %in% colnames(clinical_meta)) {
        stat_ctrl <- get_cont_stat(var, "Control")
        stat_hipe <- get_cont_stat(var, "HIPE")
        stat_lipe <- get_cont_stat(var, "LIPE")
        stat_plus <- get_cont_stat(var, "PLUS")

        # Test across groups (Kruskal-Wallis nonparametric ANOVA)
        val_all <- as.numeric(clinical_meta[[var]])
        valid_idx <- !is.na(val_all)

        if (sum(valid_idx) > 4) {
            kw_test <- kruskal.test(val_all[valid_idx] ~ factor(Arm)[valid_idx])
            p_val <- kw_test$p.value
        } else {
            p_val <- NA
        }

        clinical_res_group <- rbind(clinical_res_group, list(
            Variable = var,
            Control = paste0(stat_ctrl, " [Mean (SD)]"),
            HIPE = paste0(stat_hipe, " [Mean (SD)]"),
            LIPE = paste0(stat_lipe, " [Mean (SD)]"),
            PLUS = paste0(stat_plus, " [Mean (SD)]"),
            P_Value = p_val
        ))
    }
}

for (var in clinical_vars_cat) {
    if (var %in% colnames(clinical_meta)) {
        val_all <- clinical_meta[[var]]
        valid_idx <- !is.na(val_all)

        if (sum(valid_idx) > 0) {
            tbl <- table(val_all[valid_idx], Arm[valid_idx])

            if (nrow(tbl) > 1 && ncol(tbl) > 1) {
                suppressWarnings({
                    f_test <- fisher.test(tbl, simulate.p.value = TRUE, B = 2000)
                    p_val <- f_test$p.value
                })
            } else {
                p_val <- NA
            }

            lvl_to_report <- rownames(tbl)[length(rownames(tbl))]

            get_cat_stat <- function(group) {
                count <- sum(val_all[Arm == group] == lvl_to_report, na.rm = TRUE)
                total <- sum(!is.na(val_all[Arm == group]))
                perc <- ifelse(total > 0, count / total * 100, 0)
                sprintf("%d (%.1f%%)", count, perc)
            }

            clinical_res_group <- rbind(clinical_res_group, list(
                Variable = paste0(var, " (=", lvl_to_report, ")"),
                Control = paste0(get_cat_stat("Control"), " [n (%)]"),
                HIPE = paste0(get_cat_stat("HIPE"), " [n (%)]"),
                LIPE = paste0(get_cat_stat("LIPE"), " [n (%)]"),
                PLUS = paste0(get_cat_stat("PLUS"), " [n (%)]"),
                P_Value = p_val
            ))
        }
    }
}

clinical_res_group$P_Value <- round(clinical_res_group$P_Value, 4)
write.csv(clinical_res_group, "Results/Clinical_Baseline_Group_Differences_sex.csv", row.names = FALSE)


# --- Lipidomics Baseline Differences Among Groups ---
# Using an ANOVA-like Limma F-test to assess if ANY group significantly differs from the others
design_group <- model.matrix(~ factor(Arm))
fit_group <- lmFit(exprs_base, design_group)
fit_group <- eBayes(fit_group, robust = TRUE)

# Test coefficients 2, 3, 4 (HIPE, LIPE, PLUS vs Control intercept). An overall F-test.
res_group_f <- topTable(fit_group, coef = 2:ncol(design_group), number = Inf, sort.by = "F")
res_group_f$Lipid <- rownames(res_group_f)
res_group_f$p_value <- res_group_f$P.Value
res_group_f$adj_p_value <- res_group_f$adj.P.Val

res_group_f <- merge(res_group_f, KEGG_annot, by.x = "Lipid", by.y = "VARIABLE OWL ID")

write.csv(res_group_f, "Results/DEA_Baseline_Group_Differences_sex.csv", row.names = FALSE)

sig_baseline_lipids <- sum(res_group_f$p_value < 0.05, na.rm = TRUE)
cat(sprintf("Baseline Group Differences analysis successfully completed. Found %d nominally significant baseline lipids across groups.\n", sig_baseline_lipids))
