This project is used to fit a09 data and produce all plots in fh comparison paper.

In order to use database:
1 run "pip install espressodb"
2 please download "fh_db-db.sqlite" file from [https://drive.google.com/file/d/1JpksWCjo0sZDMfOSm7ojQKNv5VDNWHXi/view?usp=sharing], and put it at same path as this txt file. 
3 change the "NAME" in "db-config.yaml" file to the correct path of "fh_db-db.sqlite" file

Data files: "a09m310_e_gA_srcs0-15.h5", "FK_Fpi_data.h5" should be put at the same path as this txt file.

Total numb of files check: 21 (including "fh_db-db.sqlite", two data files and all hidden items)

USAGE:
After put "fh_db-db.sqlite" and two data files at right place: 
If you want plots of data or fit results on data, run "data_fit_plot.py".
If you want plots about late tsep 23 fits, run "late_tsep_23.py".
If you want stability plots or spectrum plots, run "stability_plot.py" or "spectrum_plot.py"(spectrum takes some time, about 2 min). <These two files depend on database>
If you want the rest stability plots or more spectrum plots, you need to use "fit_with_database.py" to produce the fit results you need for plots and store them in the database (just make "save=True" in the fit code), then turn to "stability_plot.py" or "spectrum_plot.py", add extra blocks.

FILES:

"fit_with_database.py": This is fit code combined with EspressoDB, you can save fit results in the database named "fh_db". Two ways of fit paras input are allowed: "scatter" and "continue"
	1 "scattered" is used to fit with scattered tsep and tau list, but it cannot be marked by tmin/ tmax in database, so you need "id_num" to mark your fits.

	2 "continuous" is used to fit with tsep list within a successive section like [tmin, tmax], so you can change tmin/ tmax in a for loop to do a series fits at a time. The stability plots with varying tmin/tmax/nstates depend on this way of fit, fixing all other paras and search for all saved tmin/tmax/nstates values.

"data_fit_plot.py": This code is used to plot processed data and fit results on data for best fits of three conbinations of fit (23s, 2s, 23). <This code dose not depend on database>

"late_tsep_23.py": This code is used to plot figures related to late tsep ([10, 12, 14]) 2pt+3pt fit, including (stability plots with varying nstates and tau cut) and (stability plots with varying 2pt tmin). <This code does not depend on database>

“stability_plot.py”: Now it can print stability plots varying [2pt tmin/ 2pt tmax/ 3pt tmin/ 3pt tmax/ sum tmin/ 3pt tau cut] stability plots of 23s fit.

"spectrum_plot.py": It can print stability plots of dEn with varying 2pt tmin, and spectrum plot of 23s best fit.

FOLDERS:

"module": This folder contains all the class files that shared by all codes, including "prior_setting.py" (prior funcs of different models), "prepare_data.py" (Class Prepare_data, used to process data file), "plot.py" and "fit.py" (Class Fit, all fit funcs).

"new_plots" This is used to store all figures.

About Database "fh_db"
This is an almost empty database (with some useful fits so that you can try "stability_plot.py" and "spectrum_plot.py"), you can use "fit_with_database.py" to fill in with the fits you need for plots.

