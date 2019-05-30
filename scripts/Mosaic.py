import numpy
import pandas
from squarify import squarify
from itertools import cycle

from bokeh.plotting import figure
from bokeh.models import (Range1d, HoverTool, ColumnDataSource, CustomJS, ColorBar, LinearColorMapper,
						  NumeralTickFormatter, FixedTicker)
from bokeh.models.widgets import Select
from bokeh.colors import RGB
from bokeh.palettes import viridis
from bokeh.io.export import export_png
from bokeh.layouts import column

from .Gene_Colors import vgene_colors, vfamily_colors, jgene_colors, isotype_colors

def Mosaic_Plot(clone_df, png = None, title = "", top_clones = 5000, count_col = "Clustered", vgene_col = "VGene",
				jgene_col = "JGene", isotype_col = "Isotype", vshm_col = "V_SHM", jshm_col = "J_SHM",
				vgene_colors = vgene_colors, vfamily_colors = vfamily_colors, jgene_colors = jgene_colors,
				isotype_colors = isotype_colors, line_width = 0.3, figsize = (600, 600), hover_tooltip = True):
	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		#"sizing_mode": "scale_both",
		"x_range": Range1d(-0.1, 1.1, bounds = (-1.0, 2.0)),
		"y_range": Range1d(-0.1, 1.1, bounds = (-1.0, 2.0)),
		#"outline_line_alpha": 0.0,
		"title": title,
		"tools": "pan, wheel_zoom, box_zoom, save, reset, help",
		"active_scroll": "wheel_zoom",
		"toolbar_location": "right"
	}

	plot = figure(**figure_params)
	plot.grid.visible = False
	plot.axis.visible = False

	hover_tooltips = [("Clone ID", "@CloneID")]

	info_cols = [count_col]
	if vgene_col is not None:
		info_cols.append(vgene_col)
		hover_tooltips.append(("V Gene", "@" + vgene_col))
	if jgene_col is not None:
		info_cols.append(jgene_col)
		hover_tooltips.append(("J Gene", "@" + jgene_col))
	if isotype_col is not None:
		info_cols.append(isotype_col)
		hover_tooltips.append(("Isotype", "@" + isotype_col))
	if vshm_col is not None:
		info_cols.append(vshm_col)
		hover_tooltips.append(("V Gene SHM", "@" + vshm_col + "{(0.00%)}"))
	if jshm_col is not None:
		info_cols.append(jshm_col)
		hover_tooltips.append(("J Gene SHM", "@" + jshm_col + "{(0.00%)}"))

	if hover_tooltip:
		hover_tool = HoverTool(point_policy = "snap_to_data", tooltips = hover_tooltips)
		plot.add_tools(hover_tool)

	mosaic_df = clone_df[info_cols]
	mosaic_df = mosaic_df.sort_values([count_col], ascending = [False])

	if top_clones:
		mosaic_df = mosaic_df.head(top_clones)

	total_area = float(mosaic_df[count_col].sum())
	mosaic_df["Clone_Frequencies"] = mosaic_df[count_col].astype(float) / total_area

	hover_tooltips.append(("Clone Frequency", "@Clone_Frequencies{(0.00%)}"))

	mosaic_rects = squarify(mosaic_df["Clone_Frequencies"].tolist(), 0.0, 0.0, 1.0, 1.0)
	#Add half width/height to x/y position for center points
	mosaic_df["x"] = [rect["x"] + rect["dx"] / 2.0 for rect in mosaic_rects]
	mosaic_df["y"] = [rect["y"] + rect["dy"] / 2.0 for rect in mosaic_rects]
	mosaic_df["width"] = [rect["dx"] for rect in mosaic_rects]
	mosaic_df["height"] = [rect["dy"] for rect in mosaic_rects]

	#By default there is no legend text, since colors are alternating and non-informative
	mosaic_df["legend"] = ""
	mosaic_df["Empty_Legend"] = ""

	alternating_colors = [RGB(102, 194, 165), RGB(252, 141, 98), RGB(141, 160, 203)]
	alt2_color_cycle = cycle(alternating_colors[0:2])
	alt3_color_cycle = cycle(alternating_colors)
	mosaic_df["alternating2_colors"] = [next(alt2_color_cycle) for _ in mosaic_rects]
	mosaic_df["alternating3_colors"] = [next(alt3_color_cycle) for _ in mosaic_rects]
	#Default color scheme is alternating 3 colors
	mosaic_df["fill_color"] = mosaic_df["alternating3_colors"]

	#Set up various mosaic coloring options and associated legends
	color_select_options = ["Alternating (2)", "Alternating (3)"]
	if vgene_col in mosaic_df.columns:
		mosaic_df["vgene_colors"] = mosaic_df[vgene_col].map(vgene_colors)
		vfamilies = mosaic_df[vgene_col].str.split("-").str[0]
		mosaic_df["vfamily_colors"] = vfamilies.map(vfamily_colors)
		color_select_options.append("V Gene")
		color_select_options.append("V Family")
		mosaic_df["VGene_Legend"] = mosaic_df[vgene_col]
		mosaic_df["VFamily_Legend"] = mosaic_df[vgene_col].str.split("-").str[0]
	if jgene_col in mosaic_df.columns:
		mosaic_df["jgene_colors"] = mosaic_df[jgene_col].map(jgene_colors)
		color_select_options.append("J Gene")
		mosaic_df["JGene_Legend"] = mosaic_df[jgene_col]
	if isotype_col in mosaic_df.columns:
		mosaic_df["isotype_colors"] = mosaic_df[isotype_col].map(isotype_colors)
		color_select_options.append("Isotype")
		mosaic_df["Isotype_Legend"] = mosaic_df[isotype_col]

	#Using viridis as a quantitative heatmap color scheme for SHM values
	#The SHM values are binned into 180 groups; viridis in >180 bins uses some values twice, which pandas.cut can't use
	shm_viridis = list(viridis(180))
	colorbar_tick_formatter = NumeralTickFormatter(format = "0.00%")

	if vshm_col in mosaic_df.columns:
		vshm_min = mosaic_df[vshm_col].min()
		vshm_max = mosaic_df[vshm_col].max()
		#Use pandas.cut to bin the V gene SHM values into the heatmap colors
		mosaic_df["vshm_colors"] = pandas.cut(mosaic_df[vshm_col], bins = 180, labels = shm_viridis)
		color_select_options.append("V Gene SHM")

		vshm_color_mapper = LinearColorMapper(palette = shm_viridis, low = vshm_min, high = vshm_max)
		vshm_ticks = FixedTicker(ticks = numpy.linspace(vshm_min, vshm_max, 8))
		vshm_colorbar = ColorBar(color_mapper = vshm_color_mapper, location = (0, 0), name = "vshm_colorbar",
								 label_standoff = 12, formatter = colorbar_tick_formatter, ticker = vshm_ticks)
		plot.add_layout(vshm_colorbar, "right")

	if jshm_col in mosaic_df.columns:
		jshm_min = mosaic_df[jshm_col].min()
		jshm_max = mosaic_df[jshm_col].max()
		#Use pandas.cut to bin the J gene SHM values into the heatmap colors
		mosaic_df["jshm_colors"] = pandas.cut(mosaic_df[jshm_col], bins = 180, labels = shm_viridis)
		color_select_options.append("J Gene SHM")

		jshm_color_mapper = LinearColorMapper(palette = shm_viridis, low = jshm_min, high = jshm_max)
		jshm_ticks = FixedTicker(ticks = numpy.linspace(jshm_min, jshm_max, 8))
		jshm_colorbar = ColorBar(color_mapper = jshm_color_mapper, location = (0, 0), name = "jshm_colorbar",
								 label_standoff = 12, formatter = colorbar_tick_formatter, ticker = jshm_ticks)
		plot.add_layout(jshm_colorbar, "right")

	mosaic_source = ColumnDataSource(mosaic_df)

	plot.rect(x = "x", y = "y", width = "width", height = "height", fill_color = "fill_color", legend = "legend",
			  line_color = "black", line_width = line_width, source = mosaic_source)

	#By default, the plot legend and ColorBar should be turned off (since the color is repeating and uninformative)
	plot.legend[0].visible = False
	vshm_colorbar = plot.select("vshm_colorbar")[0]
	jshm_colorbar = plot.select("jshm_colorbar")[0]
	vshm_colorbar.visible = False
	jshm_colorbar.visible = False

	if png is not None:
		export_png(plot, png)

	change_args = {
		"source": mosaic_source,
		"legend_obj": plot.legend[0],
		"vshm_colorbar_obj": vshm_colorbar,
		"jshm_colorbar_obj": jshm_colorbar
	}
	change_rect_color = CustomJS(args = change_args, code = """
		var selection = cb_obj.value.toLowerCase();
		var new_color_array;
		var new_legend_array;

		if(selection.indexOf("v gene shm") !== -1) {
			new_color_array = source.data["vshm_colors"];
			new_legend_array = source.data["Empty_Legend"];
			legend_obj.visible = false;
			vshm_colorbar_obj.visible = true;
			jshm_colorbar_obj.visible = false;
		} else if(selection.indexOf("j gene shm") !== -1) {
			new_color_array = source.data["jshm_colors"];
			new_legend_array = source.data["Empty_Legend"];
			legend_obj.visible = false;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = true;
		} else if(selection.indexOf("v gene") !== -1) {
			new_color_array = source.data["vgene_colors"];
			new_legend_array = source.data["VGene_Legend"];
			legend_obj.visible = true;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		} else if(selection.indexOf("v family") !== -1) {
			new_color_array = source.data["vfamily_colors"];
			new_legend_array = source.data["VFamily_Legend"];
			legend_obj.visible = true;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		} else if(selection.indexOf("j gene") !== -1) {
			new_color_array = source.data["jgene_colors"];
			new_legend_array = source.data["JGene_Legend"];
			legend_obj.visible = true;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		} else if(selection.indexOf("isotype") !== -1) {
			new_color_array = source.data["isotype_colors"];
			new_legend_array = source.data["Isotype_Legend"];
			legend_obj.visible = true;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		} else if(selection.indexOf("2") !== -1) {
			new_color_array = source.data["alternating2_colors"];
			new_legend_array = source.data["Empty_Legend"];
			legend_obj.visible = false;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		} else {
			new_color_array = source.data["alternating3_colors"];
			new_legend_array = source.data["Empty_Legend"];
			legend_obj.visible = false;
			vshm_colorbar_obj.visible = false;
			jshm_colorbar_obj.visible = false;
		}

		var fill_color = source.data["fill_color"];
		var legend = source.data["legend"];
		for(idx = 0; idx < fill_color.length; idx++) {
			fill_color[idx] = new_color_array[idx];
			legend[idx] = new_legend_array[idx];
		}
		source.change.emit();
	""")

	patch_coloring_select = Select(title = "Color by:", options = color_select_options, value = "Alternating (3)",
						   callback = change_rect_color)

	plot_layout = column(patch_coloring_select, plot)

	return plot_layout
