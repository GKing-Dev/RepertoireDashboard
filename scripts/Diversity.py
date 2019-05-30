import numpy

from bokeh.plotting import figure
from bokeh.models import Range1d, BasicTickFormatter
from bokeh.colors import RGB
from bokeh.io.export import export_png

def Shannon_Wiener_Index(clone_counts):
	"""Calculates the Shannon-Wiener index of diversity given an iterable of clone frequencies.

	Parameters
	----------
	clone_counts: iterable of floats
		Frequencies of all clones/members of a population

	Returns
	----------
	sw_index: float
		The Shannon-Wiener index value
	"""

	total_clone_counts = float(sum(clone_counts))
	sw_index = 0.0

	for clone in clone_counts:
		clone_freq = clone / total_clone_counts
		sw_index -= (clone_freq * numpy.log(clone_freq))

	return sw_index

def Hill_Diversity_Index(clone_counts, N = (0.0, 10.0), step = 0.1):
	"""Calculates the Hill Diversity index/indices; will return the indices from orders 0 to 10 by default.

	Parameters
	----------
	clone_counts: iterable of floats
		Frequencies of all clones/members of a population
	N: tuple of (float, float)
		Start and stop point for orders to calculate the diversity index from (both inclusive); default is (0.0, 10.0)
	step: float
		Step value to create the range from the start and stop point in N; default is 0.1

	Returns
	----------
	hill_index/hill_indices: float or list of floats
		The Hill Diversity index/indices for the given repertoire
	"""

	if hasattr(N, "__iter__"):
		if len(N) != 2 or N[1] <= N[0]:
			raise ValueError("If N is an iterable it must be of length 2 for the start/end orders.")
		chunks = numpy.floor(N[1] / step) + 1
		orders = numpy.linspace(start = N[0], stop = N[1], num = chunks)
	else:
		orders = [N]

	total_clone_counts = float(sum(clone_counts))
	hill_indices = []

	for order in orders:
		hill_index = 0.0
		if order == 0.0:
			#Hill index at zero is simply the species richness (total number of clones)
			hill_index = len(clone_counts)
		elif order == 1.0:
			#Hill index at one is the exponential of the Shannon-Wiener index
			sw_index = Shannon_Wiener_Index(clone_counts)
			hill_index = numpy.exp(sw_index)
		else:
			for clone in clone_counts:
				clone_freq = clone / total_clone_counts
				hill_index += (clone_freq ** order)
			order_exponent = 1.0 / (1.0 - order)
			hill_index = hill_index ** order_exponent

		order_index = (order, hill_index)
		hill_indices.append(order_index)

	if len(hill_indices) == 1:
		return hill_indices[0]
	else:
		return hill_indices

def Diversity_Plot(clone_df, png = None, title = "", count_col = "Clustered", split_col = None, line_width = 3,
				   add_control_diversities = True, figsize = (1000, 700)):
	"""Creates a plot comparing clonal repertoire diversity rates, using the Hill Diversity metric.

	Parameters
	----------
	clone_df: pandas DataFrame
		DataFrame of the repertoire(s) to plot
	png: str
		Title of the output PNG filename or None if none should be made; default is None
	title: str
		Title of the output graph; default is ""
	count_col: str
		Column name in clone_df of the clone counts/frequencies; default is "Clustered"
	split_col: str
		Column separating various repertoire subsets in clone_df or None if single repertoire; default is None
	line_width: int
		Width for the plot lines
	add_control_diversities: bool
		Whether to add lines for control diversities of artificial polarity; default is True
	figsize: tuple of (int, int)
		The width and height of the output plot; default is (100, 700)

	Returns
	----------
	plot: bokeh figure
		The figure object for the diversity plot
	"""

	figure_params = {
		"plot_width": figsize[0],
		"plot_height": figsize[1],
		"x_range": Range1d(0, 10),
		"y_axis_type": "log",
		"title": title,
		"tools": "save, help",
		"toolbar_location": "right"
	}

	plot = figure(**figure_params)
	plot.xgrid.grid_line_alpha = 0.0
	plot.xaxis.axis_label = "Order (N)"
	plot.yaxis.axis_label = "Hill Diversity Constant"
	plot.yaxis.formatter = BasicTickFormatter()

	#If comparing multiple samples, add the sample column to split on to the DataFrame
	if split_col is not None:
		diversity_df = clone_df[[count_col, split_col]]

		samples = []
		diversity_dfs = []
		for sample, df in diversity_df.groupby([split_col]):
			samples.append(sample)
			diversity_dfs.append(df)

	else:
		samples = ["Repertoire"]
		diversity_dfs = [clone_df[[count_col]]]

	sample_colors = (RGB(30, 160, 120), RGB(220, 90, 0), RGB(120, 110, 180), RGB(230, 40, 140))
	for sample, df, line_color in zip(samples, diversity_dfs, sample_colors[:len(samples)]):
		hill_indices = Hill_Diversity_Index(df[count_col])
		n_orders = [i[0] for i in hill_indices]
		order_diversities = [i[1] for i in hill_indices]

		#ADD MORE LINE STYLES (dotted, etc.)
		plot.line(x = n_orders, y = order_diversities, color = line_color, line_width = line_width, legend = sample)

	if add_control_diversities:
		total_clones = max([len(i) for i in diversity_dfs])
		total_counts = max([df[count_col].sum() for df in diversity_dfs])

		#Very highly polarized data creates a sample in which the top 20 clones are 20% of the total by prevalence
		top20_20_data = [total_counts * 0.2 / 20] * 20
		top20_20_data += [total_counts * 0.8 / (total_clones - 20) for _ in range(total_clones - 20)]
		#Highly polarized data has the top 20 clones at 15% of the total
		top20_15_data = [total_counts * 0.15 / 20] * 20
		top20_15_data += [total_counts * 0.85 / (total_clones - 20) for _ in range(total_clones - 20)]
		#Moderately polarized data has the top 20 clones at 10% of the total
		top20_10_data = [total_counts * 0.1 / 20] * 20
		top20_10_data += [total_counts * 0.9 / (total_clones - 20) for _ in range(total_clones - 20)]
		#Lowly polarized data has the top 20 clones at 5% of the total
		top20_5_data = [total_counts * 0.05 / 20] * 20
		top20_5_data += [total_counts * 0.95 / (total_clones - 20) for _ in range(total_clones - 20)]

		top20_20_diversities = [i[1] for i in Hill_Diversity_Index(top20_20_data)]
		top20_15_diversities = [i[1] for i in Hill_Diversity_Index(top20_15_data)]
		top20_10_diversities = [i[1] for i in Hill_Diversity_Index(top20_10_data)]
		top20_5_diversities = [i[1] for i in Hill_Diversity_Index(top20_5_data)]
		plot.line(x = n_orders, y = top20_20_diversities, color = RGB(160, 200, 230), alpha = 0.8, line_dash = (12,),
				  line_width = line_width, legend = "Very Highly Polarized (Top 20 Clones 20%)")
		plot.line(x = n_orders, y = top20_15_diversities, color = RGB(30, 120, 180), alpha = 0.8, line_dash = (12,),
				  line_width = line_width, legend = "Highly Polarized (Top 20 Clones 15%)")
		plot.line(x = n_orders, y = top20_10_diversities, color = RGB(180, 220, 140), alpha = 0.8, line_dash = (12,),
				  line_width = line_width, legend = "Moderately Polarized (Top 20 Clones 10%)")
		plot.line(x = n_orders, y = top20_5_diversities, color = RGB(50, 160, 40), alpha = 0.8, line_dash = (12,),
				  line_width = line_width, legend = "Lowly Polarized (Top 20 Clones 5%)")

	if png is not None:
		export_png(plot, png)

	return plot
