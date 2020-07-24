These are some help notes for you to use Formfactor_fit code.

There are three files(Formfactor_fit.ipynb, Spectrum_plots.ipynb, Read_me.txt) and one folder(module) in the same path as this txt file.



Formfactor_fit.ipynb

The main function is in Formfactor_fit.ipynb, you can run this notebook in jupyter environment. You can change parameters of fit and some other choices(whether chained or whether set a start point p0 for lsqfit) in the main function. 

If you want to fit with different data file, except for changing the 'file_name' and 'file_path' in the main function, you also need to change some parameters in the 'Class Prepare_data' according to the structure of data dictionary.

There is the 'plot' part in Formfactor_fit.ipynb, you can plot the data for estimation or plot fit results on data. But before doing plot, you need to run the main function to do a fit first.

If you need to use Djangle database, you need to put the project folder('my_project') in the same path as this txt file, and data file should be put inside the project folder('my_project').

If you do not use Djangle database, you need to delete the line "os.chdir(os.getcwd() + '/my_project')" in the 'Set up' part of both Formfactor_fit.ipynb and Spectrum_plots.ipynb, and data file should be put in the same path as this txt file.



Spectrum_plots.ipynb

This file is used to make energy spectrum plots, in 'Theoretical prediction of energy spectrum' part, you need to set proton mass and box size L according to which data file you plot.

'Fit results' part will read the fit results of En from database and check whether they are stable with varying parameter. 

You need to set values of parameters in both 'Fit results' and 'Prior and spectrum plots' parts to select the fit results you want.



module folder

Inside this folder, there are two .py file(prior_setting.py and plot.py). 

In prior_setting.py there is a function named 'prior', you can change prior of fit here.

In plot.py there is a class named 'Plot', you can change the parameters and style of plots that you get after running Formfactor_fit.ipynb in plot.py.