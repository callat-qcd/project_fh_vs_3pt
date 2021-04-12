# READ ME

This project is used to fit the three-point vector and axia-vector correlation functions on the a09m310 data.  It also produces most plots in the paper, [arXiv:2104.xxxxx](https://arxiv.org/abs/2104.xxxx)

## Preparation

### Data files
The correlation functions can be obtained with
```
wget https://a51.lbl.gov/~callat/published_results/a09m310_e_gA_srcs0-15_all_tau.h5
```
The data file includes the current insertion time in between the source and sink, as well as the "out-of-time" data.

After obtaining the file, either rename it, or create a softlink in the root directory of this repository:
```
ln -s a09m310_e_gA_srcs0-15_all_tau.h5 a09m310_e_gA_srcs0-15.h5
```

### Using the database:
An sqlite database of fit results for this project can be obtained
```
wget https://a51.lbl.gov/~callat/published_results/ga_excited_states_a09m310_db.sqlite
```
which can be used to generate plots.  Also, one can create a new database and perform fits.  

The present code is tested with (the existing database is known to not work on newer versions of gvar - this will be updated at some point in the future)
- [gvar](https://github.com/gplepage/gvar), version 11.5.2
- [lsqfit](), version 11.5.3


To make use of the database, you need to install [Espressodb](https://arxiv.org/abs/1912.03580).

```
pip install espressodb
```
Edit the "NAME" in the `db-config.yaml` file to match the downloaded sqlite file (including a relative or full path you create).  Then, to set up the db
```
python manage.py makemigrations
python manage.py migrate
```


## Figures in the paper

With the downloaded sqlite file, one can make these plots in the paper with the corresponding code:
- Fig. 1 and 2: `python fit_on_data_plot.py`
- Fig. 3: `python plot_eg_results.py`

For other plots, the plotting functions are contained in the `stability_plot.py` file.  You can run
```
python stability_plot.py -h
```
to see how to generate the various plots.

- Fig. 4 t_min/t_max: `python stability_plot.py`
- Fig. 4 t_even: `python stability_plot.py --tmin_max_stability --t_even`
- Fig. 4 t_even: `python stability_plot.py --tmin_max_stability --t_odd
- Fig. 5 excited state contamination: `python stability_plot.py`


for those plots that made by "stability_plot.py", inside the .py file, there will be a specific function to make each plot, the corresponding function is in the () below, if you would like to make these plots, you should put the corresponding functions into the main function of "stability_plot.py" and run it

all plots will be saved in "new_plots" folder

```
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
