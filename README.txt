This project is used to fit a09 data and produce all plots in fh comparison paper.

Numb of files check: 

Data files: "a09m310_e_gA_srcs0-15.h5", "FK_Fpi_data.h5" should be put at the same path as this txt file.

FILES: (13 files)

"fit_with_database.py": This is fit code combined with EspressoDB, you can save fit results in the database named "fh_db". Two ways of fit paras input are allowed: "scatter" and "continue"
	"scattered" is used to fit with scattered tsep and tau list, but it cannot be marked by tmin/ tmax in database, so you 		need "id_num" to mark your fits.

	"continuous" is used to fit with tsep list within a successive section like [tmin, tmax], so you can change tmin/ tmax in 		a for loop to do a series fits at a time. The stability plots with varying tmin/tmax/nstates depend on this way of fit, 			fixing all other paras and search for all saved tmin/tmax/nstates values.

"data_fit_plot.py": This code is used to plot processed data and fit results on data for best fits of three conbinations of fit (23s, 2s, 23). <This code dose not depend on database>

"late_tsep_23.py": This code is used to plot figures related to late tsep ([10, 12, 14]) 2pt+3pt fit, including (stability plots with varying nstates and tau cut) and (stability plots with varying 2pt tmin). <This code does not depend on database>

“stability_plot.py”: This code is used to produce all stability plots.

FOLDERS: (4 folders)

"module": This folder contains all the class files that shared by all codes, including "prior_setting.py" (prior funcs of different models), "prepare_data.py" (Class Prepare_data, used to process data file), "plot.py" and "fit.py" (Class Fit, all fit funcs).

"new_plots" This is used to store all figures. <IF THERE IS NO "new_plots" folder, please create one first>

About Database "fh_db"
This is an almost empty database, you can use "fit_with_database.py" to fill in with the fits you need for plots.

