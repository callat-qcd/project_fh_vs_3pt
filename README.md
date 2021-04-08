# READ ME

This project is used to fit a09 data and produce all plots in fh comparison paper.

## Prepare

download a09 h5 file from somewhere and save it in the same folder as this README:

```
a09m310_e_gA_srcs0-15.h5
```

In order to use database:

```
1. run "pip install espressodb"

2. download "fh_db-db.sqlite" file and save it in the same folder as this README (this db file contains all fits that you needed to make every plot)

wget <file_from_andre_desktop>

3. change the "NAME" in "db-config.yaml" file to the correct path of "fh_db-db.sqlite" file 

4. run "python manage.py makemigrations"

5. run "python manage.py migrate"
```

## Figures in the paper


you can make these plots in the paper with corresponding code

for those plots that made by "fit_on_data_plot.py" and "excited_state.py", just need to run the .py file

for those plots that made by "stability_plot.py", inside the .py file, there will be a specific function to make each plot, the corresponding function is in the () below, if you would like to make these plots, you should put the corresponding functions into the main function of "stability_plot.py" and run it

all plots will be saved in "new_plots" folder

```
Fig.1: fit_on_data_plot.py 

Fig.3: stability_plot.py ( tmin_tmax_combined() , both_even_tmax() and both_odd_tmax() )

Fig.4: excited_state.py 

Fig.5: fit_on_data_plot.py

Fig.6: stability_plot.py ( late_tsep_23_tau_inc() )

Fig.7: stability_plot.py ( late_tsep_23_E0() )

Fig.8: stability_plot.py ( ga_summary() )

Fig.10: fit_on_data_plot.py 

Fig.11: stability_plot.py ( pt2_tmin_2s() , sum_tmin_2s() )

Fig.12: stability_plot.py ( pt2_tmax_2s() )

Fig.13: stability_plot.py ( pt2_tmin_23() , pt3_tmin_23() )

Fig.15: stability_plot.py ( tau_cut_23() )

Fig.16: stability_plot.py ( pt2_tmin_23s() , pt3_tmin_23s() and sum_tmin_23s() )

Fig.17: stability_plot.py ( tau_cut_23s() )

Fig.18: stability_plot.py ( prior_width_23s() )

```

## Do the fit with database

```
fit_with_database.py

you can change parameters to do any fit you want, and choose to store results in the database or not.
```

## About things in module folder

```
fit.py : all fit functions and process about fits

p0.py : best_p0 to speed up fit process

plot.py : all plot functions

prepare_data.py : deal with the h5 data file, do average and output a dict

prior_setting.py : all prior models

```