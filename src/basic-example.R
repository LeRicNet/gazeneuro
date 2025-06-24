library(gazeneuro)
library(tidyverse)

nifti_data <- preload_nifti_data("../2025_ATPC_Study/data/img/1.3.12.2.1107.5.2.30.26451.2007113015563474695781235.0.0.0.nii.gz")
gaze_data <- read_csv('../2025_ATPC_Study/data/gaze_data/records-drzmexhret-dab29f3dcfc6a92938d61db1d396e492c4590042625149988a148b9a95ccddf1.csv')
z_axis <- read_csv('../2025_ATPC_Study/data/scroll_index_atpc.csv')

# need to filter z-axis to match gaze data session. all gaze_data files are saves as 'records-<clientID>-<tracking_session>'.
# note: z_axis$narrative_session == gaze_data$tracking_session
z_axis <- z_axis %>%
  filter(narrative_session == gaze_data$tracking_session[1])

# 3. Integrate the data
integrated <- integrate_all_gaze_points(gaze_data, z_axis)

# 4. Visualize
plot_slice_with_all_gaze(nifti_data, integrated, slice_num = 13)
