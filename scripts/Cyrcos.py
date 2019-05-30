import numpy
from itertools import combinations

from bokeh.plotting import figure
from bokeh.models import Range1d, ColumnDataSource
from bokeh.palettes import all_palettes
from bokeh.io import save, show
from bokeh.embed import components

class Cyrcos_Repertoire_Comparison_Plot(object):
	def __init__(self, clone_dfs, title = "", top_clones = None, normalize_segments = True, gap_size = 10,
				 start_pos = "top", clockwise = True, offset_segments = None, segment_face_colors = "Category10",
				 segment_outline_colors = None, fade_segments = True, clone_col = "CloneID", count_col = "Clustered",
				 sample_col = "Sample", figsize = (1000, 1000)):
		"""Creates a Circos-like Chord graph for comparing multiple immune repertoire clonotype profiles."""

		#Plot visual aspect definitions
		segment_width = 0.07 #Thickness of the circle arc segments
		offset_shift_amount = 0.1 #Increase in radius for segments to offset
		min_segment_alpha = 0.2

		#Gather the samples and their respective DataFrames
		df_cols = [clone_col, count_col]
		self.samples = []
		comparison_dfs = []
		if isinstance(clone_dfs, dict): #Input is a dictionary of {sample_name: DataFrame}
			for sample in clone_dfs:
				self.samples.append(sample)
				clone_df = clone_dfs[sample][df_cols].sort_values([count_col], ascending = [False])

				if top_clones is not None:
					clone_df = clone_df.head(top_clones)

				comparison_dfs.append(clone_df)

		else: #Input is a single DataFrame with all samples; sample_col gives the sample to split on
			df_cols.append(sample_col)
			for sample, df in clone_dfs[df_cols].groupby([sample_col]):
				self.samples.append(sample)
				clone_df = df.sort_values([count_col], ascending = [False])

				if top_clones is not None:
					clone_df = clone_df.head(top_clones)

				comparison_dfs.append(clone_df)

		self.total_samples = len(self.samples)

		#Create the figure plot
		self.Create_Plot(title = title, figsize = figsize)

		sample_clone_counts = [len(df) for df in comparison_dfs]
		total_gap_size = gap_size * self.total_samples
		total_segment_len = 360 - total_gap_size

		if normalize_segments:
			self.segment_lengths = [total_segment_len / self.total_samples] * self.total_samples
		else:
			total_clones = sum(sample_clone_counts)
			self.segment_lengths = [clones / total_clones * total_segment_len for clones in sample_clone_counts]

		#Convert description of the first segment's start location to the angular position in degrees:
		start_pos = start_pos.lower()
		if "top" in start_pos or "north" in start_pos:
			start_position = 90
		elif "right" in start_pos or "east" in start_pos:
			start_position = 0
		elif "bottom" in start_pos or "south" in start_pos:
			start_position = -90
		elif "left" in start_pos or "west" in start_pos:
			start_position = -180
		else:
			start_position = 90

		#Shift the start position to align the first gap's center with the start position.
		self.direction = "clock" if clockwise else "anticlock"
		start_position -= gap_size / 2 if clockwise else 0

		radius = 0.8
		self.inner_radii = [radius] * self.total_samples

		if offset_segments is not None:
			if isinstance(offset_segments, int):
				self.inner_radii[offset_segments] += offset_shift_amount

			elif hasattr(offset_segments, "__iter__") and not isinstance(offset_segments, str):
				for seg in offset_segments:
					self.inner_radii[seg] += offset_shift_amount

		self.outer_radii = [r + segment_width for r in self.inner_radii]

		if isinstance(segment_face_colors, str):
			if segment_face_colors in all_palettes:
				self.segment_face_colors = all_palettes[segment_face_colors][self.total_samples]
			else:
				self.segment_face_colors = [segment_face_colors] * self.total_samples
		elif hasattr(segment_face_colors, "__iter__") and not isinstance(segment_face_colors, str):
			if len(segment_face_colors) == self.total_samples:
				self.segment_face_colors = segment_face_colors
			else:
				raise IndexError("List provided to segment_face_colors is not the same length as total segments!")
		else:
			print("Warning: segment_face_colors should be a colormap or list of colors for each segment!")
			print("Defaulting to \"Category10\"...")
			self.segment_face_colors = all_palettes["Category10"][self.total_samples]

		if segment_outline_colors is None:
			self.segment_outline_colors = ["transparent"] * self.total_samples
		elif isinstance(segment_outline_colors, str):
			self.segment_outline_colors = [segment_outline_colors] * self.total_samples
		elif hasattr(segment_outline_colors, "__iter__") and not isinstance(segment_outline_colors, str):
			self.segment_outline_colors = segment_outline_colors
		else:
			raise TypeError("segment_outline_colors should be a color name, list of colors, or None!")

		#Add the repertoire circle segments to the figure
		self.Create_Segments(start_position, gap_size, clockwise, fade_segments, min_alpha = min_segment_alpha)

		#Add rank and segment position columns to the clone DataFrames
		for idx, _ in enumerate(comparison_dfs):
			#Rank clones by total count / frequency (largest clone being rank 0, second largest rank 1, etc.)
			clone_ranks = comparison_dfs[idx][count_col].rank(method = "first", ascending = False).astype(int)
			comparison_dfs[idx]["Rank"] = clone_ranks - 1

			#Convert the rank to a relative position from 0.0 to 1.0 for placement along the segments
			comparison_dfs[idx]["Position"] = comparison_dfs[idx]["Rank"] / len(comparison_dfs[idx])

		#Iterate through all combinations of two samples and create the links from clone to clone
		for comb in combinations(range(self.total_samples), 2):
			idx1 = comb[0]
			idx2 = comb[1]
			sample1 = self.samples[idx1]
			sample2 = self.samples[idx2]
			sample_df1 = comparison_dfs[idx1].set_index([clone_col])
			sample_df2 = comparison_dfs[idx2].set_index([clone_col])

			#Join the samples into a DataFrame containing only the shared clones and their positions
			joined_df = sample_df1.join(sample_df2, how = "inner", lsuffix = "1", rsuffix = "2")

			#Calculate the angular position for the clones (segment start location + clone position * segment length)
			#If the plot is drawn clockwise, subtract the segment start instead of adding it
			seg1_start = -self.segment_starts[idx1] if self.direction == "clock" else self.segment_starts[idx1]
			seg2_start = -self.segment_starts[idx2] if self.direction == "clock" else self.segment_starts[idx2]
			pos1 = joined_df["Position1"] * self.segment_lengths[idx1] + seg1_start
			pos2 = joined_df["Position2"] * self.segment_lengths[idx2] + seg2_start

			#Convert the positions to the start and end xy coordinates
			inner_rad1 = self.inner_radii[idx1]
			inner_rad2 = self.inner_radii[idx2]
			xy1 = self.Angle_to_XY(angles = pos1, radius = inner_rad1)
			xy2 = self.Angle_to_XY(angles = pos2, radius = inner_rad2)
			#Convert the list of (x, y) tuples to a list of xs and a list of ys
			xs1, ys1 = zip(*xy1)
			xs2, ys2 = zip(*xy2)

			link_data = {
				"x0": xs1,
				"y0": ys1,
				"x1": xs2,
				"y1": ys2
			}
			link_source = ColumnDataSource(link_data)

			#Plot the links matching clone positions between repertoires; control points are the center of the circle
			self.plot.quadratic(x0 = "x0", y0 = "y0", x1 = "x1", y1 = "y1", cx = 0.5, cy = 0.5, source = link_source,
								color = "black", line_width = 1)

	def Create_Plot(self, title, figsize):
		plot_params = {
			"plot_width": figsize[0],
			"plot_height": figsize[1],
			"x_range": Range1d(-0.5, 1.5, bounds = (-1.0, 2.0)),
			"y_range": Range1d(-0.5, 1.5, bounds = (-1.0, 2.0)),
			"title": title,
			"tools": "pan, wheel_zoom, box_zoom, save, reset, help",
			"active_scroll": "wheel_zoom",
			"toolbar_location": "right",
			"outline_line_alpha": 0.0
		}
		self.plot = figure(**plot_params)
		self.plot.grid.visible = False
		self.plot.axis.visible = False

	def Create_Segments(self, start_position, gap_size, clockwise, fade_segments, fade_steps = 1000, min_alpha = 0.01):
		self.segment_starts = []
		self.segment_ends = []

		if clockwise:
			segment_deltas = [-seg_len for seg_len in self.segment_lengths]
			gap_delta = -gap_size
		else:
			segment_deltas = [seg_len for seg_len in self.segment_lengths]
			gap_delta = gap_size

		cur_position = start_position

		segment_border_starts = []
		segment_border_ends = []

		for seg_len in segment_deltas:
			segment_border_starts.append(cur_position)
			cur_position += seg_len
			segment_border_ends.append(cur_position)
			cur_position += gap_delta

		border_source_dict = {
			"start_angle": segment_border_starts,
			"end_angle": segment_border_ends,
			"inner_radius": self.inner_radii,
			"outer_radius": self.outer_radii,
			"line_color": self.segment_outline_colors
		}
		border_data = ColumnDataSource(data = border_source_dict)

		self.borders = self.plot.annular_wedge(x = 0.5, y = 0.5, start_angle = "start_angle", end_angle = "end_angle",
											   inner_radius = "inner_radius", outer_radius = "outer_radius",
											   fill_color = None, line_color = "line_color", line_width = 1,
											   direction = self.direction, start_angle_units = "deg",
											   end_angle_units = "deg", source = border_data)

		cur_position = start_position
		cur_seg_starts = []
		cur_seg_ends = []

		if fade_segments:
			self.segment_alphas = []

			#Extend the radius and fill color lists to account for the extra alpha segments:
			inner_seg_radii = [r for r in self.inner_radii for _ in range(fade_steps)]
			outer_seg_radii = [r for r in self.outer_radii for _ in range(fade_steps)]
			self.segment_face_colors = [seg_color for seg_color in self.segment_face_colors for _ in range(fade_steps)]

			for segment_delta in segment_deltas:
				alpha_segment_delta = segment_delta / fade_steps
				alpha_delta = (1.0 - min_alpha) / fade_steps
				cur_alpha = 1.0

				self.segment_starts.append(cur_position)

				for _ in range(fade_steps):
					cur_seg_starts.append(cur_position)
					cur_position += alpha_segment_delta
					cur_seg_ends.append(cur_position)

					self.segment_alphas.append(cur_alpha)
					cur_alpha -= alpha_delta

				self.segment_ends.append(cur_position)
				cur_position += gap_delta

			cur_legend = [label for label in self.samples for _ in range(fade_steps)]

		else:
			self.segment_alphas = [1.0] * self.total_samples
			inner_seg_radii = self.inner_radii
			outer_seg_radii = self.outer_radii

			for segment_delta in segment_deltas:
				self.segment_starts.append(cur_position)
				cur_seg_starts.append(cur_position)

				cur_position += segment_delta
				self.segment_ends.append(cur_position)
				cur_seg_ends.append(cur_position)

				cur_position += gap_delta

			cur_legend = self.samples

		source_data_dict = {
			"start_angle": cur_seg_starts,
			"end_angle": cur_seg_ends,
			"inner_radius": inner_seg_radii,
			"outer_radius": outer_seg_radii,
			"fill_color": self.segment_face_colors,
			"fill_alpha": self.segment_alphas,
			"legend": cur_legend
		}
		source_data = ColumnDataSource(data = source_data_dict)

		self.wedges = self.plot.annular_wedge(x = 0.5, y = 0.5, direction = self.direction, start_angle = "start_angle",
											  end_angle = "end_angle", inner_radius = "inner_radius",
											  outer_radius = "outer_radius", fill_color = "fill_color",
											  fill_alpha = "fill_alpha", line_color = None, start_angle_units = "deg",
											  end_angle_units = "deg", legend = "legend", source = source_data)

	def Angle_to_XY(self, angles, radius, angles_in_degrees = True, offset = (0.5, 0.5)):
		"""Converts an angular position to X, Y coordinates.

		Parameters
		----------
		angles: int/float or iterable of int/float
			Angle(s) to convert to X, Y coordinate(s).
		radius: int/float
			Radius of the circle for which the X, Y coordinates map to.
		angles_in_degrees: bool
			Whether the provided angle(s) are in degrees or radians; default is True.
		offset: tuple of (float, float)
			The x and y location of the center of the circle; default is (0.5, 0.5).

		Returns
		----------
		xy_coords: tuple of (float, float) or list of tuples (float, float)
			X and Y coordinate(s) of the angles provided on a circle with provided radius.
		"""

		if angles_in_degrees:
			x = radius * numpy.sin(numpy.deg2rad(angles)) + offset[0]
			y = radius * numpy.cos(numpy.deg2rad(angles)) + offset[1]
		else:
			x = radius * numpy.sin(angles) + offset[0]
			y = radius * numpy.cos(angles) + offset[1]

		#Return a tuple of (x, y) if a single angle was provided, or a list of (x, y) tuples if multiple angles
		if hasattr(x, "__iter__") and hasattr(y, "__iter__"):
			xy_coords = list(zip(x, y))
		else:
			xy_coords = (x, y)

		return xy_coords

	def Show(self):
		"""Call Bokeh show function to display the current plot in a browser window."""

		show(self.plot)

	def Get_Plot_Components(self):
		"""Call Bokeh components function to get the HTML script and div elements representing the current plot.

		Returns
		----------
		plot_script: str
			The HTML code for the Bokeh JavaScript plot to display.
		plot_div: str
			The HTML code for the div element used to place the plot in a page.
		"""

		plot_script, plot_div = components(self.plot)
		return plot_script, plot_div

	def Save_HTML(self, filename, title = "Repertoire Comparison Cyrcos Graph"):
		"""Calls Bokeh save function to save an HTML file containing the current plot.

		Parameters
		----------
		filename: str
			The desired filename / filepath for the output HTML file.
		title: str
			A title to use for the HTML file; default is "Repertoire Comparison Cyrcos Graph".
		"""

		save(self.plot, filename = filename, title = title)
