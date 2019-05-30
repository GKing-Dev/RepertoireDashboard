import numpy
import pandas

from bokeh.io import save, show, output_file
from bokeh.layouts import layout

from scripts.Diversity import Diversity_Plot
from scripts.Cyrcos import Cyrcos_Repertoire_Comparison_Plot
from scripts.UpSet import Repertoire_Upset_Plot
from scripts.Mosaic import Mosaic_Plot
from scripts.Gene_Plots import VJ_Gene_Plot, Burtin_VGene_SHM_Plot
from scripts.Clone_Stats import Violin_SHM_Plot, CDR_Length_Histogram_Plot
from scripts.Gene_Colors import vgene_colors, vfamily_colors, jgene_colors, isotype_colors

def Repertoire_Dashboard(clone_dfs, filename = None, title = "Repertoire Analysis Dashboard", plot_title_prefix = "",
						 mosaic_top_clones = 5000, cyrcos_top_clones = 1000, upset_highlighted_sets = None,
						 clone_col = "CloneID", vgene_col = "VGene", jgene_col = "JGene", isotype_col = "Isotype",
						 count_col = "Clustered", vshm_col = "V_SHM", jshm_col = "J_SHM", cdr_col = "CDR3_AA",
						 sample_col = None, sizing_mode = "scale_width", show_plots = True, bokeh_resources = "cdn"):
	"""Creates an interactive dashboard HTML page displaying all the comparative visualizations.

	Parameters
	----------
	clone_dfs: pandas DataFrame or dict of {str: DataFrame}
		Input repertoire(s) with sample names; input should be formatted as a dict of sample name: DataFrame or a single
		concatenated pandas DataFrame with a column sample_col specifying the samples of origin
	filename: str or None
		Name for the saved output HTML file, or None if user wishes to save manually; default is None
	title: str
		Page title for the output HTML file; default is "Repertoire Analysis Dashboard"
	plot_title_prefix: str
		Optional prefix added to the title of each plot, for example to prefix titles with a Donor name; default is ""
	mosaic_top_clones: int
		Limit for the total clones to display for Mosaic plots (over 5000 is often visually jarring); default is 5000
	cyrcos_top_clones: int
		Limit for the total clones to display for Cyrcos plots; default is 1000
	upset_highlighted_sets: list of tuples or None
		Specific shared sample sets to highlight in different colors for the UpSet comparison plot; default is None
	clone_col: str
		Header / name for the column containing the sample clone IDs; default is "CloneID"
	vgene_col: str
		Header / name for the column containing the sample V genes; default is "VGene"
	jgene_col: str
		Header / name for the column containing the sample J genes; default is "JGene"
	isotype_col: str
		Header / name for the column containing the sample isotypes; default is "Isotype"
	count_col: str
		Header / name for the column containing the sample clone counts or frequencies; default is "Clustered"
	vshm_col: str
		Header / name for the column containing the sample clone V gene SHMs; default is "V_SHM"
	jshm_col: str
		Header / name for the column containing the sample clone J gene SHMs; default is "J_SHM"
	cdr_col: str
		Header / name for the column containing the sample clone CDR3 amino acid sequences; default is "CDR3_AA"
	sample_col: str or None
		Header / name for the column containing sample names if all samples are in one DataFrame; default is "Sample"
	sizing_mode: str
		How to scale the plots in the dashboard layout (see Bokeh layout function); default is "scale_width"
	show_plots: bool
		Whether the dashboard page will be immediately shown upon creation; default is True
	bokeh_resources: str
		BokehJS resource location used for the dashboard (see Bokeh output_file documentation); default is "cdn"
		"cdn" gets the required files from the Bokeh CDN (requires internet connection)
		"inline" adds all necessary stylesheets and scripts to the HTML page itself

	Returns
	----------
	dashboard: bokeh nested layout of Column and Row
		The output dashboard Layout object representing the final plots and their placements
	"""

	#Set up output file if user wants to save the dashboard page
	if filename is not None:
		output_file(filename = filename, title = title, mode = bokeh_resources)

	repertoire_cols = [clone_col, vgene_col, jgene_col, isotype_col, count_col, vshm_col, jshm_col, cdr_col]

	if isinstance(clone_dfs, dict):
		if sample_col is None:
			sample_col = "Sample"

		dfs = []
		for sample in clone_dfs:
			clone_df = clone_dfs[sample][repertoire_cols]
			clone_df[sample_col] = sample
			dfs.append(clone_df)

		comparison_df = pandas.concat(dfs, ignore_index = True)

	elif isinstance(clone_dfs, pandas.DataFrame):
		if sample_col is not None:
			repertoire_cols.append(sample_col)
			comparison_df = clone_dfs[repertoire_cols]
		else:
			raise KeyError("No sample-name column header was found in the repertoire DataFrame!")

	#############################################
	##    Paired V-J Gene Usage Donut Plots    ##
	#############################################
	vj_gene_plots = []
	for sample, df in comparison_df.groupby([sample_col]):
		plot_title = "{0} {1} Paired V-J Gene Usage".format(plot_title_prefix, sample)
		vj_gene_plot = VJ_Gene_Plot(df, title = plot_title, vgene_col = vgene_col, jgene_col = jgene_col,
									count_col = count_col, vgene_colors = vgene_colors, jgene_colors = jgene_colors,
									vfamily_colors = vfamily_colors)
		vj_gene_plots.append(vj_gene_plot)

	#############################################
	##        V/J Gene SHMs Violin Plot        ##
	#############################################
	vj_shm_plot = Violin_SHM_Plot(comparison_df, title = plot_title_prefix + " Gene SHM Levels", vshm_col = vshm_col,
								  jshm_col = jshm_col, split_col = sample_col)

	#############################################
	## Repertoire Clone Frequency Mosaic Plots ##
	#############################################
	mosaic_plots = []
	for sample, df in comparison_df.groupby([sample_col]):
		plot_title = "{0} {1} Clonotype Frequencies Mosaic".format(plot_title_prefix, sample)
		mosaic_plot = Mosaic_Plot(df, title = plot_title, top_clones = mosaic_top_clones, vgene_col = vgene_col,
								  jgene_col = jgene_col, isotype_col = isotype_col, count_col = count_col,
								  vshm_col = vshm_col, jshm_col = jshm_col, vgene_colors = vgene_colors,
								  jgene_colors = jgene_colors, vfamily_colors = vfamily_colors,
								  isotype_colors = isotype_colors)
		mosaic_plots.append(mosaic_plot)

	#############################################
	##      Clonal V Gene SHM Burtin Plot      ##
	#############################################
	clonal_vgene_shm_plot = Burtin_VGene_SHM_Plot(comparison_df, title = plot_title_prefix + " Clonal V Gene Mean SHM",
												  vgene_col = vgene_col, vshm_col = vshm_col, split_col = sample_col,
												  vfamily_colors = vfamily_colors)

	#############################################
	##        Repertoire Diversity Plot        ##
	#############################################
	diversity_plot = Diversity_Plot(comparison_df, title = plot_title_prefix + " Repertoire Diversity & Polarization",
									count_col = count_col, split_col = sample_col)

	#############################################
	##  CDR3 Amino Acid Length Histogram Plot  ##
	#############################################
	cdr_len_plot = CDR_Length_Histogram_Plot(comparison_df, title = plot_title_prefix + " CDR3 Length Spectratype",
											 cdr_col = cdr_col, split_col = sample_col)

	#############################################
	## Shared Repertoire Clonotypes UpSet Plot ##
	#############################################
	upset_plot = Repertoire_Upset_Plot(comparison_df, title = plot_title_prefix + " Shared Clone Set UpSet Plot",
									   clone_col = clone_col, sample_col = sample_col,
									   highlighted_sets = upset_highlighted_sets)

	#############################################
	## Shared Clone Rank/Frequency Circos Plot ##
	#############################################
	cyrcos_plot = Cyrcos_Repertoire_Comparison_Plot(comparison_df, title = " Shared Repertoire Clonal Frequency",
													top_clones = cyrcos_top_clones, clone_col = clone_col,
													count_col = count_col, sample_col = sample_col)

	dashboard_layout = [[upset_plot.plots_grid], [vj_shm_plot, clonal_vgene_shm_plot], [cdr_len_plot, diversity_plot]]

	#Arrange the mosaic and V-J gene plots into two columns
	mosaic_plots = [list(plots) for plots in numpy.array_split(mosaic_plots, numpy.ceil(len(mosaic_plots) / 2))]
	vj_gene_plots = [list(plots) for plots in numpy.array_split(vj_gene_plots, numpy.ceil(len(vj_gene_plots) / 2))]

	dashboard_layout += mosaic_plots
	dashboard_layout += vj_gene_plots

	dashboard_layout += [[cyrcos_plot.plot]]
	dashboard = layout(children = dashboard_layout, sizing_mode = sizing_mode)

	if show_plots:
		show(dashboard)
	else:
		save(dashboard)

if __name__ == "__main__":
	df = pandas.read_csv("data/Donor_Clones.txt", sep = "\t", usecols = ["CloneID", "Clustered", "VGene", "JGene", "Isotype", "V_SHM", "J_SHM", "CDR3_AA", "Sample"],
						 dtype = {"CloneID": int, "Clustered": int, "VGene": str, "JGene": str, "Isotype": str, "V_SHM": float, "J_SHM": float, "CDR3_AA": str, "Sample": str})
	Repertoire_Dashboard(df, filename = "Repertoire_Dashboard.html", plot_title_prefix = "Donor", sample_col = "Sample")
