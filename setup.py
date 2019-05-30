from setuptools import setup

setup(
	name = "RepertoireDashboard",
	version = "1.0",
	author = "Greg King",
	author_email = "greg.king@utexas.edu",
	description = "Automated immune repertoire comparison dashboard",
	long_description = "",
	url = "https://github.com/GKing-Dev/RepertoireDashboard",
	scripts = [
		"scripts/Repertoire_Dashboard.py",
		"scripts/Diversity",
		"scripts/Cyrcos",
		"scripts/UpSet",
		"scripts/Mosaic",
		"scripts/Gene_Plots",
		"scripts/Clone_Stats",
		"scripts/Gene_Color"
	],
	install_requires = ["pandas", "numpy", "bokeh", "scipy", "squarify"],
	packages = ["RepertoireDashboard"],
	license = "MIT",
	keywords = ["immune repertoire", "antibody", "repertoire", "dashboard", "visualization", "immunology"]
)
