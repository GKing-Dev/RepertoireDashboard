import numpy
import json
import math

from bokeh.plotting import figure
from bokeh.models import Range1d, HoverTool, ColumnDataSource, CustomJS
from bokeh.models.widgets import Select
from bokeh.colors import RGB
from bokeh.io import save
from bokeh.io.export import export_png
from bokeh.layouts import column, layout, Spacer

from .Gene_Colors import vgene_colors, vfamily_colors, jgene_colors

def VJ_Gene_Plot(clone_df, png = None, title = "", vgene_col = "VGene", jgene_col = "JGene", count_col = "Clustered",
				 vgene_colors = vgene_colors, vfamily_colors = vfamily_colors, jgene_colors = jgene_colors,
				 vj_gap = 0.008, vgene_gap = 0.0, line_width = 0.4, figsize = (800, 800), hover_tooltip = True):
	"""Creates a donut (??) chart for prevalence of all V/J gene pairs in a Repertoire.

	Parameters
	----------
	clone_df: pandas DataFrame

	Returns
	----------
	script: str
	div: str
	"""

	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		#"sizing_mode": "scale_both",
		"x_range": Range1d(-0.5, 1.5, bounds = (-1.5, 2.5)),
		"y_range": Range1d(-0.5, 1.5, bounds = (-1.5, 2.5)),
		#"outline_line_alpha": 0.0,
		"title": title,
		"tools": "pan, wheel_zoom, box_zoom, tap, save, reset, help",
		"active_scroll": "wheel_zoom",
		"toolbar_location": "right"
	}

	plot = figure(**figure_params)
	plot.grid.visible = False
	plot.axis.visible = False

	if hover_tooltip:
		hover_tool = HoverTool(tooltips = [("Gene", "@legend"), ("Percent", "@percent{(0.00%)}")],
									   point_policy = "snap_to_data")
		plot.add_tools(hover_tool)

	gene_df = clone_df[[vgene_col, jgene_col, count_col]].groupby([vgene_col, jgene_col]).agg({count_col: sum})
	#Sort by V gene ascending, then J gene ascending
	gene_df = gene_df.sort_index()
	gene_df = gene_df.reset_index()

	total_vgenes = len(gene_df[vgene_col].drop_duplicates())
	total_gapsize = total_vgenes * vgene_gap
	remaining_size = 360.0 - float(total_gapsize)
	gap_size = float(vgene_gap)

	total_counts = gene_df[count_col].sum()
	gene_df["Arc_Length"] = gene_df[count_col] / total_counts * remaining_size
	#Starting at 90 degrees (top center of the circle) plus half the gap size
	#cur_v_start = -90.0 + (gap_size / 2.0)
	cur_v_start = 90.0 + (gap_size / 2.0)

	v_start_angles = []
	v_end_angles = []
	vgene_facecolors = []
	vgene_hover_colors = []
	vfamily_facecolors = []
	vfamily_hover_colors = []
	v_legend_text = []
	v_legend_percent = []

	j_start_angles = []
	j_end_angles = []
	jgene_facecolors = []
	jgene_hover_colors = []
	j_legend_text = []
	j_legend_percent = []

	for vgene in gene_df[vgene_col].drop_duplicates():
		cur_vgene_df = gene_df[gene_df[vgene_col] == vgene]
		vfamily = vgene.split("-")[0]

		vgene_color = vgene_colors[vgene]
		vgene_hover_color = vgene_color.darken(0.05)
		vfamily_color = vfamily_colors[vfamily]
		vfamily_hover_color = vfamily_color.darken(0.05)

		v_arc_length = cur_vgene_df["Arc_Length"].sum()
		cur_v_end = cur_v_start + v_arc_length

		v_start_angles.append(cur_v_start)
		v_end_angles.append(cur_v_end)

		vgene_facecolors.append(vgene_color)
		vgene_hover_colors.append(vgene_hover_color)
		vfamily_facecolors.append(vfamily_color)
		vfamily_hover_colors.append(vfamily_hover_color)

		v_legend_text.append(vgene)
		cur_vgene_counts = cur_vgene_df["Clustered"].sum()
		v_legend_percent.append(cur_vgene_counts / total_counts)

		cur_j_start = cur_v_start
		for jgene, jgene_arc_length in zip(cur_vgene_df[jgene_col], cur_vgene_df["Arc_Length"]):
			cur_j_end = cur_j_start + jgene_arc_length

			jgene_color = jgene_colors[jgene]
			jgene_hover_color = jgene_color.darken(0.05)

			j_start_angles.append(cur_j_start)
			j_end_angles.append(cur_j_end)

			jgene_facecolors.append(jgene_color)
			jgene_hover_colors.append(jgene_hover_color)

			cur_j_start = cur_j_end

			j_legend_text.append(jgene)
			cur_jgene_counts = cur_vgene_df[cur_vgene_df[jgene_col] == jgene]["Clustered"].sum()
			j_legend_percent.append(cur_jgene_counts / cur_vgene_counts)

		cur_v_start = cur_v_end + gap_size

	v_wedge_data = {
		"start_angle": v_start_angles,
		"end_angle": v_end_angles,
		"fill_color": vgene_facecolors,
		"legend": v_legend_text,
		"percent": v_legend_percent,
		"vgene_facecolors": vgene_facecolors,
		"vfamily_facecolors": vfamily_facecolors,
		"hover_fill_color": vgene_hover_colors,
		"vgene_hover_colors": vgene_hover_colors,
		"vfamily_hover_colors": vfamily_hover_colors
	}
	v_source = ColumnDataSource(v_wedge_data)

	v_inner_rad = 0.4
	v_outer_rad = 0.692

	plot.annular_wedge(x = 0.5, y = 0.5, start_angle = "start_angle", end_angle = "end_angle",
					   fill_color = "fill_color", selection_fill_color = "fill_color",
					   nonselection_fill_color = "fill_color", selection_fill_alpha = 1.0,
					   nonselection_fill_alpha = 0.2, hover_fill_color = "hover_fill_color", inner_radius = v_inner_rad,
					   outer_radius = v_outer_rad, line_color = "black", line_width = line_width, source = v_source,
					   legend = "legend", start_angle_units = "deg", end_angle_units = "deg")

	j_wedge_data = {
		"start_angle": j_start_angles,
		"end_angle": j_end_angles,
		"fill_color": jgene_facecolors,
		"legend": j_legend_text,
		"percent": j_legend_percent,
		"hover_fill_color": jgene_hover_colors
	}

	j_source = ColumnDataSource(j_wedge_data)

	j_inner_rad = v_outer_rad + vj_gap
	j_outer_rad = j_inner_rad + 0.15

	plot.annular_wedge(x = 0.5, y = 0.5, start_angle = "start_angle", end_angle = "end_angle",
					   fill_color = "fill_color", selection_fill_color = "fill_color",
					   nonselection_fill_color = "fill_color", selection_fill_alpha = 1.0,
					   nonselection_fill_alpha = 0.2, hover_fill_color = "hover_fill_color", inner_radius = j_inner_rad,
					   outer_radius = j_outer_rad, line_color = "black", line_width = line_width, source = j_source,
					   legend = "legend", start_angle_units = "deg", end_angle_units = "deg")

	if png is not None:
		export_png(plot, png)

	change_v_color = CustomJS(args = {"source": v_source}, code = """
		var selection = cb_obj.value;
		var new_color_array;
		var new_hover_array;
		if(selection.toLowerCase().indexOf("gene") !== -1) {
			new_color_array = source.data["vgene_facecolors"];
			new_hover_array = source.data["vgene_hover_colors"];
		} else {
			new_color_array = source.data["vfamily_facecolors"];
			new_hover_array = source.data["vfamily_hover_colors"];
		}
		var fill_color = source.data["fill_color"];
		var hover_fill_color = source.data["hover_fill_color"];
		for(idx = 0; idx < fill_color.length; idx++) {
			fill_color[idx] = new_color_array[idx];
			hover_fill_color[idx] = new_hover_array[idx];
		}
		source.change.emit();
	""")

	v_data_color_by = Select(title = "Color by:", options = ["V Gene", "V Family"],
							 value = "V Gene", callback = change_v_color)

	plot_layout = column(v_data_color_by, plot)
	return plot_layout

def Burtin_VGene_SHM_Plot(clone_df, png = None, title = "", vgene_col = "VGene", vshm_col = "V_SHM", split_col = None,
						  vfamily_colors = vfamily_colors, label_arc = 20, figsize = (900, 900)):
	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		"x_axis_type": None,
		"y_axis_type": None,
		"x_range": Range1d(-45, 45, bounds = (-50, 50)),
		"y_range": Range1d(-45, 45, bounds = (-50, 50)),
		"title": title,
		"tools": "pan, wheel_zoom, box_zoom, save, reset, help",
		"active_scroll": "wheel_zoom",
		"toolbar_location": "right",
		"background_fill_color": RGB(216, 216, 216)
	}

	plot = figure(**figure_params)
	plot.grid.visible = False
	plot.axis.visible = False

	label_offset = 90 #Offset the SHM % labels to the top of the plot
	plot_data_degrees = 360 - label_arc
	initial_angle = label_offset + label_arc / 2
	ending_angle = label_offset + 360 - label_arc / 2
	plot_inner_rad = 10
	plot_outer_rad = 35
	plot_thickness = plot_outer_rad - plot_inner_rad

	df_cols = [vgene_col, vshm_col]
	#If comparing multiple samples, add the sample column to split on to the DataFrame
	if split_col is not None:
		df_cols.append(split_col)

	vgene_shm_df = clone_df[df_cols].sort_values([vgene_col])

	total_vgenes = len(vgene_shm_df[vgene_col].drop_duplicates())
	vgene_arc_degrees = plot_data_degrees / total_vgenes

	#Create and color arc backgrounds by V family
	vgene_family_df = vgene_shm_df[[vgene_col]].drop_duplicates().reset_index(drop = True)
	vgene_family_df["VFamily"] = vgene_family_df[vgene_col].str.split("-").str[0]
	vgene_family_df["fill_color"] = vgene_family_df["VFamily"].map(vfamily_colors)
	vfamily_arc_length = plot_data_degrees / total_vgenes
	vgene_family_df["start_angle"] = vgene_family_df.index * vfamily_arc_length + initial_angle
	vgene_family_df["end_angle"] = vgene_family_df["start_angle"] + vfamily_arc_length

	vfamily_source = ColumnDataSource(vgene_family_df)
	plot.annular_wedge(x = 0, y = 0, start_angle = "start_angle", end_angle = "end_angle", fill_color = "fill_color",
					   inner_radius = plot_inner_rad, outer_radius = plot_outer_rad, line_color = None,
					   source = vfamily_source, start_angle_units = "deg", end_angle_units = "deg")

	if split_col in vgene_shm_df:
		vgene_shm_dfs = [sample_df_tup for sample_df_tup in vgene_shm_df.groupby([split_col])]

		samples = []
		grouped_vgene_shm_dfs = []
		for sample, df in vgene_shm_dfs:
			samples.append(sample)
			grouped_vgene_shm_dfs.append(df.groupby([vgene_col])[vshm_col].agg({"mean"}))

		#Add the V genes that may be present in one sample but not in the current one
		all_vgenes = vgene_shm_df[vgene_col].drop_duplicates().tolist()
		grouped_vgene_shm_dfs = [df.reindex(all_vgenes).reset_index() for df in grouped_vgene_shm_dfs]

		vshm_min = min([df["mean"].min() for df in grouped_vgene_shm_dfs])
		vshm_max = max([df["mean"].max() for df in grouped_vgene_shm_dfs])

	else:
		grouped_vgene_shm_df = vgene_shm_df.groupby([vgene_col])[vshm_col].agg({"mean"}).reset_index()
		grouped_vgene_shm_df = grouped_vgene_shm_df.sort_values([vgene_col]).reset_index(drop = True)

		vshm_min = grouped_vgene_shm_df["mean"].min()
		vshm_max = grouped_vgene_shm_df["mean"].max()

		samples = ["All"]
		grouped_vgene_shm_dfs = [grouped_vgene_shm_df]

	#Create the labels and radial axis lines for the SHM data
	shm_labels = ["{0:.1%}".format(shm) for shm in numpy.linspace(vshm_min, vshm_max, 7)]
	shm_label_radii = numpy.linspace(plot_inner_rad, plot_outer_rad, 7)
	plot.circle(x = 0, y = 0, radius = shm_label_radii, fill_color = None, line_color = "white")
	plot.text(x = 0, y = shm_label_radii[1:], text = shm_labels[1:], text_font_size = "10pt",
			  text_align = "center", text_baseline = "middle")

	#Create line-width annular wedges to separate V genes
	sep_angles = numpy.linspace(initial_angle, ending_angle, total_vgenes + 1)
	sep_inner_radius = plot_inner_rad - 1
	sep_outer_radius = plot_outer_rad + 1
	plot.annular_wedge(x = 0, y = 0, start_angle = sep_angles, end_angle = sep_angles, fill_color = None,
					   inner_radius = sep_inner_radius, outer_radius = sep_outer_radius, line_color = "black",
					   start_angle_units = "deg", end_angle_units = "deg")

	#Gene text labels; text angle location is the midpoint of the V gene separation lines
	text_radius = plot_outer_rad + 3.5
	text_radian_locs = numpy.deg2rad((sep_angles[1:] + sep_angles[:-1]) / 2)
	text_x = text_radius * numpy.cos(text_radian_locs)
	text_y = text_radius * numpy.sin(text_radian_locs)
	#Angle the text based on the position around the circle; reverse the left half so the text isn't upside-down
	mid_graph_radian = numpy.deg2rad(label_offset + 180)
	text_angles = [rad if rad > mid_graph_radian else rad + numpy.pi for rad in text_radian_locs]
	plot.text(x = text_x, y = text_y, text = vgene_family_df[vgene_col], angle = text_angles,
			  text_font_size = "10pt", text_align = "center", text_baseline = "middle")

	#Finally draw the bars and legend for the mean SHM values for all clones of a specific V gene
	total_samples = len(grouped_vgene_shm_dfs)
	vgene_arc_radians = numpy.deg2rad(vgene_arc_degrees)
	bar_width = vgene_arc_radians / (total_samples + 1)
	spacer_width = bar_width / (total_samples + 1)
	sample_colors = (RGB(60, 60, 60), RGB(130, 40, 40), RGB(60, 60, 130), RGB(10, 50, 100), RGB(150, 100, 20))
	sample_label_ys = numpy.linspace(-total_samples, total_samples, total_samples)
	arc_starts = text_radian_locs - (vgene_arc_radians / 2) + spacer_width

	for sample, cur_df in enumerate(grouped_vgene_shm_dfs):
		bar_start_angles = arc_starts + sample * (bar_width + spacer_width)
		bar_end_angles = bar_start_angles + bar_width
		cur_df["Normalized_SHM"] = cur_df["mean"] / vshm_max

		shm_bars = cur_df["Normalized_SHM"] * plot_thickness + plot_inner_rad
		plot.annular_wedge(x = 0, y = 0, start_angle = bar_start_angles, end_angle = bar_end_angles, line_color = None,
					   inner_radius = plot_inner_rad, outer_radius = shm_bars, fill_color = sample_colors[sample])

		if total_samples > 1:
			plot.rect(x = -2, y = sample_label_ys[sample], width = 2.5, height = 1.5, color = sample_colors[sample])
			plot.text(x = 0, y = sample_label_ys[sample], text = {"value": samples[sample]}, text_font_size = "10pt",
					  text_baseline = "middle")

	if png is not None:
		export_png(plot, png)

	return plot
