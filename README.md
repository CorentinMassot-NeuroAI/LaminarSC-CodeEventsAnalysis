# LaminarSC-CodeEventsAnalysis
MATLAB code for estimating and displaying events in spiking activity from laminar data recorded in the Superior Colliculus.


‘analysis_#’ functions
they are the main display and analysis functions of the data across all sessions
they call subfunctions or rely on saved results of subfunctions named ‘compute_#’.

analysis_vmi: displays VMI

analysis_peak_latencyfrompeak: displays onset of burst activity during visual epoch

analysis_delay_avg: displays activity during delay epoch

analysis_onset_buildup_avg: displays onsets of buildup and burst events during pre-saccadic epoch


‘compute_#’ functions: compute features of the spiking activity for each session
usually save the results in a separate matrix for each session (in directory results/).

‘get_#’ functions: dedicated subfunctions called by ‘analysis_#’ or ‘compute_#’ functions

get_inflection_2pwlr: computes inflection point using 2 piece-wise linear fit between 2 points.

‘plot_#’ functions: functions to plot results 

plot_stats_depths_v: plot results with their statistical estimation across depth and for a vertical display


Reference: please cite paper:  Massot C., Jagadisan U.K. & Gandhi N.J., Sensorimotor transformation elicits systematic patterns of activity along the dorsoventral extent of the superior colliculus in the macaque monkey, (to appear in) Communications Biology, 2019.

