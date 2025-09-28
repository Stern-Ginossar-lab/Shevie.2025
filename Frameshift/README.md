# Frameshift Simulation

The evaluation of potential frameshift signals in Ribo-Seq data is carried out with **four main scripts**:

---

## 1. **NO_frameshift_simulation.Rmd**  
**Purpose:**  
Estimate the expected fraction of reads in each ribosomal frame (`F`, `Fp1`, `Fm1`) under a **no-frameshift null model**.  

**Approach:**  
- Resample codon-level counts across a range of window sizes.  
- Report quantiles, means, and standard deviations for each frame and for the ratio `Fp1/F`.

---

## 2. **calc_fractions_on_sites.Rmd**  
**Purpose:**  
Scan coding sequences (CDS) for candidate **frameshift-associated windows**.  

**Approach:**  
- Identify sites such as:  
  - Slippery motifs (e.g., `UUUC`).  
  - Codon pairs encoding the same amino acid in 0 and +1 frames (“SS”).  
  - User-defined sets of positions.  
- Summarize frame-specific read counts within each window.  
- Score each window against the no-frameshift null expectations.

---

## 3. **frameshift_simulation.R**  
**Purpose:**  
Model the effect of introducing different **percentages of frameshifted reads** into Ribo-Seq–like data.  

**Approach:**  
- Mix synthetic frameshifted reads with normal reads at varying percentages.  
- Recalculate frame distributions and Z-scores.  
- Evaluate how signals deviate from the no-frameshift simulation and from candidate window expectations.

---

## 4. **plot_simulation_and_actual_fractions.Rmd**  
**Purpose:**  
Compare **observed window-level frameshift scores** from experimental data to **simulated distributions** generated at controlled frameshift percentages.  

**Approach:**  
- Plot simulated distributions across a range of frameshift levels.  
- Overlay actual observed values to assess which simulated percentage best matches real data.  
- Export comparative figures (density plots, PDFs).
