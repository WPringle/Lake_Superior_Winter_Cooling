## This directory contains the code to create the figures from the paper
scripts are labeled with the figure numbers from the main paper and supplementary that they generate. 

some scripts use different combinations of input arguments to generate the desired figure. For example, `fig4+5+S6+S7_plot_LES_TURB_eastern_southern_mooring_comparison_T_profiles_combined.jl` uses different combinations of days (15 30 or 45 60) and forcing (`coare_wind` or `direct`) arguments, e.g.:
`julia plot_..._combined.jl coare_wind 15 30` will produce fig4.

`figSX_plot_LES_TURB_eastern_southern_mooring_comparison_T_profiles.jl` will produce plots similar to `fig4+5+S6+S7_plot_LES_TURB_eastern_southern_mooring_comparison_T_profiles_combined.jl` except for a single day (2 rows of panels), rather than two days combined (4 rows of panels).
