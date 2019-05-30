import json
import pandas

from bokeh.plotting import figure
from bokeh.models import Range1d, ColumnDataSource
from bokeh.io import save, show
from bokeh.embed import components
from bokeh.layouts import gridplot, Spacer

class Repertoire_Upset_Plot(object):
	def __init__(self, clone_dfs, title = "", min_shared = 2, max_shared = None, overlap_bounds = None,
				 clone_col = "CloneID", sample_col = None, highlighted_sets = None, figsize = (1200, 900)):
		"""Creates a Repertoire comparison UpSet overlap plot."""

		df_cols = [clone_col]

		if isinstance(clone_dfs, dict):
			df_dict = {sample: clone_dfs[sample][df_cols] for sample in clone_dfs}

			comparison_df = pandas.concat(list(df_dict.values()), ignore_index = True)
			samples = [i for i in df_dict]

		else:
			df_cols.append(sample_col)

			comparison_df = clone_dfs[df_cols]
			samples = clone_dfs[sample_col].drop_duplicates().tolist()

		#Collapse by clone IDs and concatenate the multiple repertoire sample names containing each clone using JSON
		serializer = lambda x: json.dumps(x.tolist())
		grouped_df = comparison_df.groupby([clone_col])[sample_col].agg({serializer, "count"})
		grouped_df = grouped_df.rename(columns = {"<lambda>": "Samples", "count": "Sample_Count"})

		if min_shared is not None:
			grouped_df = grouped_df[grouped_df["Sample_Count"] >= min_shared]
		if max_shared is not None:
			grouped_df = grouped_df[grouped_df["Sample_Count"] <= max_shared]

		#Calculate the total number of clones shared by each combination of samples
		overlap_counts = grouped_df["Samples"].value_counts()
		if overlap_bounds is not None:
			overlap_counts = overlap_counts[overlap_counts >= overlap_bounds[0]]
			overlap_counts = overlap_counts[overlap_counts <= overlap_bounds[1]]

		total_sets = len(overlap_counts)
		total_samples = len(samples)
		MAIN_BAR_WIDTH = 0.5
		set_colors = ("#82C882", "#BEB4D2", "#FABE82", "#FFFF96", "#326EB4", "#F00082", "#BE5A14", "#646464")

		main_plot_params = {
			"plot_width": int(figsize[0] * 0.6),
			"plot_height": int(figsize[1] * 0.75),
			"x_range": Range1d(-MAIN_BAR_WIDTH, total_sets - 1 + MAIN_BAR_WIDTH),
			"title": title,
			"tools": "save, reset, help",
			"outline_line_alpha": 0.0
		}
		self.main_plot = figure(**main_plot_params)
		self.main_plot.grid.visible = False
		self.main_plot.xaxis.visible = False
		largest_overlap = overlap_counts.max()
		self.main_plot.yaxis.bounds = (0, largest_overlap)
		self.main_plot.yaxis.axis_label_text_font_size = "12pt"

		sample_sets_xs = [i for i in range(total_sets)]
		default_bar_color = "#96AAC8"
		sample_sets_colors = []
		cur_color_idx = 0
		cur_color = default_bar_color

		if highlighted_sets is not None:
			for main_bar_set in overlap_counts.index:
				for highlighted_subset in highlighted_sets:
					cur_set = json.loads(main_bar_set)
					if set(cur_set) == set(highlighted_subset):
						cur_color = set_colors[cur_color_idx]
						cur_color_idx += 1
						break
					else:
						cur_color = default_bar_color
				sample_sets_colors.append(cur_color)
		else:
			sample_sets_colors += [default_bar_color] * len(overlap_counts)

		main_bars_data = {
			"x": sample_sets_xs,
			"top": overlap_counts.tolist(),
			"color": sample_sets_colors
		}
		main_bars_source = ColumnDataSource(main_bars_data)
		self.main_bars = self.main_plot.vbar(x = "x", top = "top", width = MAIN_BAR_WIDTH, bottom = 0,
											 color = "color", source = main_bars_source)
		self.main_plot.yaxis.axis_label = "Total Shared Clones"

		sample_sets_plot_params = {
			"plot_width": int(figsize[0] * 0.6),
			"plot_height": int(figsize[1] * 0.25),
			"x_range": self.main_plot.x_range,
			"y_range": Range1d(-MAIN_BAR_WIDTH, total_samples - 1 + MAIN_BAR_WIDTH, bounds = (0, total_samples)),
			"tools": "save, reset, help",
			"background_fill_color": "#C8C896",
			"background_fill_alpha": 0.1,
			"outline_line_alpha": 0.0
		}
		self.sample_sets_plot = figure(**sample_sets_plot_params)
		self.sample_sets_plot.grid.visible = False
		self.sample_sets_plot.xaxis.visible = False
		self.sample_sets_plot.yaxis.axis_label = " "  #Helps keep plot at the same width as the main plot dimensions
		self.sample_sets_plot.yaxis.axis_label_text_font_size = "12pt"

		#Used to pad the Y axis to align properly with the main bar graph
		self.sample_sets_plot.yaxis[0].ticker = [i for i in range(total_samples)]
		self.sample_sets_plot.yaxis.major_label_overrides = {"0": str(largest_overlap)}

		#"Draw" invisible axis line / labels / ticks to match main plot width
		self.sample_sets_plot.yaxis.axis_label_text_color = None
		self.sample_sets_plot.yaxis.axis_line_color = None
		self.sample_sets_plot.yaxis.major_tick_line_color = None
		self.sample_sets_plot.yaxis.major_label_text_color = None

		sample_set_circle_radius = MAIN_BAR_WIDTH * 0.3

		sample_sets_data = {
			"x": sample_sets_xs * total_samples,
			"y": [i for i in range(total_samples) for _ in range(total_sets)]
		}
		sample_sets_source = ColumnDataSource(sample_sets_data)
		self.sample_sets_plot.circle(x = "x", y = "y", radius = sample_set_circle_radius, fill_color = "#787878",
									 line_color = None, alpha = 0.5, source = sample_sets_source)

		#Get the total clones per sample for the total clone counts bar graph
		clone_counts = comparison_df[sample_col].value_counts()

		#Create the linked circles that mark the compared samples
		ypos_to_sample = clone_counts.reset_index()["index"].to_dict()
		sample_to_ypos = {ypos_to_sample[key]: key for key in ypos_to_sample}
		cur_set_color = "black"
		cur_color_idx = 0
		for x_pos, cur_set in enumerate(overlap_counts.index):
			sample_set = json.loads(cur_set)
			cur_set_count = len(sample_set)
			set_ys = [sample_to_ypos[sample] for sample in sample_set]

			if highlighted_sets is not None:
				for subset in highlighted_sets:
					if set(sample_set) == set(subset):
						cur_set_color = set_colors[cur_color_idx]
						cur_color_idx += 1
						break
					else:
						cur_set_color = "black"

			#Draw the bars linking the sample set circles
			min_circle_y = min(set_ys)
			max_circle_y = max(set_ys)
			self.sample_sets_plot.line(x = [x_pos, x_pos], y = [min_circle_y, max_circle_y],
									   line_width = 5, line_color = cur_set_color)
			#Draw the sample set circles
			self.sample_sets_plot.circle(x = [x_pos] * cur_set_count, y = set_ys, radius = sample_set_circle_radius,
										 fill_color = cur_set_color, line_color = None)

		largest_repertoire_clones = clone_counts.max()
		clone_bar_plot_params = {
			"plot_width": int(figsize[0] * 0.4),
			"plot_height": int(figsize[1] * 0.25),
			"x_range": Range1d(largest_repertoire_clones, 0),
			"tools": "save, reset, help",
			"outline_line_alpha": 0.0,
			"y_axis_location": "right"
		}
		self.clone_bar_plot = figure(**clone_bar_plot_params)
		self.clone_bar_plot.grid.visible = False
		self.clone_bar_plot.xaxis.axis_label_text_font_size = "12pt"
		self.clone_bar_plot.yaxis[0].ticker = [i for i in range(total_samples)]
		sample_ticks_to_labels = {tick_y: sample for tick_y, sample in enumerate(clone_counts.index.tolist())}
		self.clone_bar_plot.yaxis.major_label_overrides = sample_ticks_to_labels
		self.clone_bar_plot.yaxis.major_tick_line_color = None
		self.clone_bar_plot.yaxis.major_label_text_baseline = "middle"
		self.clone_bar_plot.yaxis.major_label_text_font_size = "12pt"

		clone_bars_data = {
			"y": [i for i in range(total_samples)],
			"right": clone_counts.tolist()
		}
		clone_bars_source = ColumnDataSource(clone_bars_data)
		self.clone_bars = self.clone_bar_plot.hbar(y = "y", right = "right", height = 0.5, left = 0,
												   color = "#96AAC8", source = clone_bars_source)
		self.clone_bar_plot.xaxis.axis_label = "Total Clones"

		top_left_spacer = Spacer(width = int(figsize[0] * 0.4), height = int(figsize[1] * 0.75))
		self.plots_grid = gridplot([top_left_spacer, self.main_plot, self.clone_bar_plot, self.sample_sets_plot],
								   ncols = 2, tools = "save, reset, help", toolbar_location = "right")

	def Show(self):
		"""Call Bokeh show function to display the current plot in a browser window."""

		show(self.plots_grid)

	def Get_Plot_Components(self):
		"""Call Bokeh components function to get the HTML script and div elements representing the current plot.

		Returns
		----------
		plot_script: str
			The HTML code for the Bokeh JavaScript plot to display.
		plot_div: str
			The HTML code for the div element used to place the plot in a page.
		"""

		plot_script, plot_div = components(self.plots_grid)
		return plot_script, plot_div

	def Save_HTML(self, filename, title = "Repertoire Comparison UpSet Graph"):
		"""Calls Bokeh save function to save an HTML file containing the current plot.

		Parameters
		----------
		filename: str
			The desired filename / filepath for the output HTML file.
		title: str
			A title to use for the HTML file; default is "Repertoire Comparison UpSet Graph".
		"""

		save(self.plots_grid, filename = filename, title = title)
