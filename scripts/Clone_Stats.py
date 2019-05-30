import numpy
import math
from scipy.stats import gaussian_kde

from bokeh.plotting import figure
from bokeh.models import Range1d, HoverTool, ColumnDataSource, NumeralTickFormatter, FixedTicker
from bokeh.io.export import export_png

def Violin_SHM_Plot(clone_df, png = None, title = "", vshm_col = "V_SHM", jshm_col = "J_SHM", split_col = None,
					quads = True, violin_width = 0.8, line_width = 0.4, figsize = (1000, 600), hover_tooltip = True):
	"""Creates a SHM violin plot that can be used to compare multiple categories in a Repertoire.

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
		"y_range": Range1d(-0.005, 0.3, bounds = (-0.01, 0.31)),
		"title": title,
		"tools": "pan, wheel_zoom, box_zoom, save, reset, help",
		"active_scroll": "wheel_zoom",
		"toolbar_location": "right"
	}

	plot = figure(**figure_params)
	plot.grid.visible = False
	plot.xaxis.minor_tick_line_color = None
	plot.xaxis.major_label_text_font_size = "10pt"
	plot.yaxis.axis_label = "V/J Gene SHM"
	plot.yaxis.major_label_text_font_size = "10pt"
	plot.yaxis.formatter = NumeralTickFormatter(format = "0.00%")

	if hover_tooltip:
		hover_tooltips = [("Mean SHM", "@mean{(0.00%)}"), ("Max SHM", "@max{(0.00%)}"),
						  ("25th Percentile", "@quantile25{(0.00%)}"), ("75th Percentile", "@quantile75{(0.00%)}")]
		hover_tool = HoverTool(point_policy = "follow_mouse", tooltips = hover_tooltips)
		plot.add_tools(hover_tool)

	shm_cols = []
	if vshm_col is not None:
		shm_cols.append(vshm_col)
	if jshm_col is not None:
		shm_cols.append(jshm_col)

	#To compare samples, add the sample column to split on to the DataFrame
	if split_col is not None:
		shm_cols.append(split_col)
		shm_df = clone_df[shm_cols]

		samples = []
		shm_dfs = []
		for sample, df in shm_df.groupby([split_col]):
			samples.append(sample)
			shm_dfs.append(df)

	else:
		samples = ["Repertoire"]
		shm_dfs = [clone_df[shm_cols]]

	vshm_violin_color = "lightgreen"
	jshm_violin_color = "slateblue"
	violin_xs = []
	violin_ys = []
	violin_colors = []
	violin_legends = []
	hover_means = []
	hover_maxes = []
	hover_25quantiles = []
	hover_75quantiles = []
	x_location_to_category = {}
	violin_x_offset = 0

	for sample, df in zip(samples, shm_dfs):
		#Create the density functions
		if vshm_col in df.columns:
			vshm_mean = df[vshm_col].mean()
			vshm_max = df[vshm_col].max()
			hover_means.append([vshm_mean])
			hover_maxes.append([vshm_max])
			hover_25quantiles.append([df[vshm_col].quantile(0.25)])
			hover_75quantiles.append([df[vshm_col].quantile(0.75)])

			y_points = numpy.linspace(0.0, vshm_max, 300)  #Create the y range of 300 points from min to max
			reversed_y_points = numpy.flipud(y_points)
			v_kernel = gaussian_kde(df[vshm_col], "scott")
			vshm_x_points = v_kernel(y_points)

			#Normalize the x range to standard width; negate V SHM points to place it on the left half of the violin
			vshm_x_points = -vshm_x_points / vshm_x_points.max() * violin_width / 2.0

			#Return to the patch starting points if a different violin is drawn for the other half, or mirror data
			if jshm_col in df.columns:
				vshm_x_points = numpy.append(vshm_x_points, abs(vshm_x_points).min())
				vshm_y_points = numpy.append(y_points, y_points.min())
			else:
				reversed_vshm_x = numpy.flipud(-vshm_x_points)
				vshm_x_points = numpy.append(vshm_x_points, reversed_vshm_x)
				vshm_y_points = numpy.append(y_points, reversed_y_points)

			violin_xs.append(vshm_x_points + violin_x_offset)
			violin_ys.append(vshm_y_points)
			violin_colors.append(vshm_violin_color)
			violin_legends.append("V Gene SHM")

		if jshm_col in df.columns:
			jshm_mean = df[jshm_col].mean()
			jshm_max = df[jshm_col].max()
			hover_means.append([jshm_mean])
			hover_maxes.append([jshm_max])
			hover_25quantiles.append([df[jshm_col].quantile(0.25)])
			hover_75quantiles.append([df[jshm_col].quantile(0.75)])

			y_points = numpy.linspace(0.0, jshm_max, 300)  #Create the y range of 300 points from min to max
			reversed_y_points = numpy.flipud(y_points)
			j_kernel = gaussian_kde(df[jshm_col], "scott")
			jshm_x_points = j_kernel(y_points)

			#Normalize the x range to standard width
			jshm_x_points = jshm_x_points / jshm_x_points.max() * violin_width / 2.0

			#Return to the patch starting points if a different violin is drawn for the other half, or mirror data
			if vshm_col in df.columns:
				jshm_x_points = numpy.append(jshm_x_points, abs(jshm_x_points).min())
				jshm_y_points = numpy.append(y_points, y_points.min())
			else:
				reversed_jshm_x = numpy.flipud(-jshm_x_points)
				jshm_x_points = numpy.append(jshm_x_points, reversed_jshm_x)
				jshm_y_points = numpy.append(y_points, reversed_y_points)

			violin_xs.append(jshm_x_points + violin_x_offset)
			violin_ys.append(jshm_y_points)
			violin_colors.append(jshm_violin_color)
			violin_legends.append("J Gene SHM")

		if quads:
			pass

		x_location_to_category[violin_x_offset] = sample
		violin_x_offset += violin_width * 1.2

	violin_data = {
		"xs": violin_xs,
		"ys": violin_ys,
		"fill_color": violin_colors,
		"legend": violin_legends,
		"mean": hover_means,
		"max": hover_maxes,
		"quantile25": hover_25quantiles,
		"quantile75": hover_75quantiles
	}
	violin_source = ColumnDataSource(violin_data)

	plot.patches(xs = "xs", ys = "ys", fill_color = "fill_color", line_color = "black", line_width = line_width,
				 legend = "legend", source = violin_source)

	#Replace / remap the X axis tickers to the categorical samples
	plot.xaxis.ticker = FixedTicker(ticks = [loc for loc in x_location_to_category])
	plot.xaxis.major_label_overrides = x_location_to_category
	plot.x_range.bounds = (min(x_location_to_category.keys()) - 1, max(x_location_to_category.keys()) + 1)

	if png is not None:
		export_png(plot, png)

	return plot

def CDR_Length_Histogram_Plot(clone_df, png = None, title = "", cdr_col = "CDR3_AA", split_col = None,
							  quantile_boundries = (0.0001, 0.9999), figsize = (800, 600)):
	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		"title": title,
		"tools": "pan, wheel_zoom, box_zoom, save, reset, help",
		"active_scroll": "wheel_zoom",
		"toolbar_location": "right"
	}

	plot = figure(**figure_params)
	plot.grid.visible = False
	plot.xaxis.minor_tick_line_color = None
	plot.xaxis.axis_label = "CDR3 Length"
	plot.xaxis.axis_label_text_font_size = "12pt"
	plot.xaxis.major_label_text_font_size = "12pt"
	plot.yaxis.axis_label = "P(x)"
	plot.yaxis.axis_label_text_font_size = "12pt"
	plot.yaxis.major_label_text_font_size = "12pt"

	#To compare samples, add the sample column to split on to the DataFrame
	if split_col is not None:
		cdr3_df = clone_df[[cdr_col, split_col]]

		samples = []
		cdr3_lens = []
		for sample, df in cdr3_df.groupby([split_col]):
			samples.append(sample)
			cdr3_lens.append(df[cdr_col].str.len())

	else:
		samples = ["Repertoire"]
		cdr3_lens = [clone_df[cdr_col].str.len()]

	bin_min = min([cdr_len_series.min() for cdr_len_series in cdr3_lens])
	bin_max = max([cdr_len_series.max() for cdr_len_series in cdr3_lens]) + 1
	bin_range = [i for i in range(bin_min, bin_max)]

	bar_colors = ["#A0C8E6", "#32A032", "#1E78B4", "#B4DC8C"]
	bar_offset = 0.0
	bar_width = 1 / len(samples)

	upper_y = 0.0

	for idx, (sample, cdr_len_series) in enumerate(zip(samples, cdr3_lens)):
		heights, lefts = numpy.histogram(cdr_len_series, density = True, bins = bin_range)

		#Ensure proper Y axis scrolling boundaries are set
		if heights.max() > upper_y:
			upper_y = heights.max()

		#Shift bars if multiple samples are being plotted
		lefts = lefts.astype(float)
		lefts += bar_offset

		bar_lefts = lefts[:-1]
		bar_rights = bar_lefts + bar_width

		plot.quad(top = heights, bottom = 0, left = bar_lefts, right = bar_rights, fill_color = bar_colors[idx],
				  line_color = None, legend = sample)

		bar_offset += bar_width

	plot.y_range.start = -0.001
	plot.y_range.end = upper_y
	plot.y_range.bounds = (-0.05, upper_y * 1.5)

	if quantile_boundries is not None:
		lower_x = clone_df[cdr_col].str.len().quantile(quantile_boundries[0])
		upper_x = clone_df[cdr_col].str.len().quantile(quantile_boundries[1])

		plot.x_range.start = lower_x
		plot.x_range.end = upper_x

	plot.x_range.bounds = (0, bin_max + 4)

	if png is not None:
		export_png(plot, png)

	return plot

def Rarefaction_Plot(align_df, png = None, title = "", cdr_col = "CDR3_AA", split_col = None, cdr_identity = 0.96,
					 steps = 50, reads = None, figsize = (800, 600), hover_tooltip = True, save_to_file = False):
	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		"title": title,
		"tools": "save, help",
		"toolbar_location": "right"
	}
	plot = figure(**figure_params)
	plot.xgrid.grid_line_color = None
	plot.xaxis.axis_label = "Total Sampled Reads"
	plot.yaxis.axis_label = "Total Clonotypes"
	plot.axis.formatter = NumeralTickFormatter(format = "0")

	tooltips = [("Total Sampled Reads", "@xs"), ("Total Clones", "@ys")]
	if hover_tooltip:
		hover_tool = HoverTool(point_policy = "snap_to_data", tooltips = tooltips, mode = "hline", names = ["rar_line"])
		plot.add_tools(hover_tool)

	#If comparing multiple samples, add the sample column to split on to the DataFrame
	if split_col is not None:
		reads_df = align_df[[cdr_col, split_col]]

		samples = []
		reads_dfs = []
		for sample, df in reads_df.groupby([split_col]):
			samples.append(sample)
			reads_dfs.append(df)

		if hover_tooltip and len(samples) > 1:
			tooltips.append(("Sample", "@sample"))

	else:
		samples = ["Repertoire"]
		reads_dfs = [align_df[[cdr_col]]]

	sample_colors = ["#1EA078", "#DC5A00", "#786EB4", "#E6288C", "#B4D28C", "#A028B4"]
	for sample, df, color in zip(samples, reads_dfs, sample_colors[:len(samples)]):
		total_reads = len(df)
		subsamp_sizes = []
		cur_total = 0

		if reads is not None:
			subsamp_steps = reads
		else:
			subsamp_steps = math.floor(total_reads / steps)

		#Create the list of all read subsample counts to clonotype
		while cur_total < total_reads:
			if cur_total != 0:
				subsamp_sizes.append(cur_total)

			cur_total += subsamp_steps

		subsamp_sizes.append(total_reads)

		subsamp_clones = []
		for n in subsamp_sizes:
			sub_read_df = df.sample(n)
			sub_total_clones = Clonotype_Usearch(sub_read_df[cdr_col], identity = cdr_identity)
			subsamp_clones.append(sub_total_clones)

		rarefaction_data = {
			"reads": subsamp_sizes,
			"clones": subsamp_clones,
			"sample": [sample if len(samples) > 1 else None] * len(subsamp_sizes)
		}
		rar_source = ColumnDataSource(rarefaction_data)

		plot.line(x = "reads", y = "clones", color = color, line_width = 3, source = rar_source,
				  legend = "sample", name = "rar_line")
		plot.scatter(x = "reads", y = "clones", color = color, source = rar_source)

		if save_to_file:
			with open(sample + "_Rarefaction_Data.txt", "w") as rarefaction_text_file:
				rarefaction_text_file.write("Reads\tClones\n")
				for read_count, clone_count in zip(subsamp_sizes, subsamp_clones):
					rarefaction_text_file.write("{0}\t{1}\n".format(read_count, clone_count))

	if png is not None:
		export_png(plot, png)

	return plot
