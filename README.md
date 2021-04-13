# READ ME

This project is used to fit the three-point vector and axia-vector correlation functions on the a09m310 data.  It also produces most plots in the paper, [arXiv:2104.05226](https://arxiv.org/abs/2104.05226)

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
- Fig. 6 excited state contamination: `python stability_plot.py`
- Fig. 7 late time 3pt fit: `python stability_plot.py --tmin_max_stability --t_large_stab`
- Fig. 8 late time 2pt sensitivity: `python stability_plot.py --tmin_max_stability --t_large_2pt`
- Fig. 9 gA summary plot: `python stability_plot.py --tmin_max_stability --ga_summary`
- Fig. 10 [in the works]:
- Fig. 11 effective mass/z0/ga plots: `python stability_plot.py`
- Fig. 12+13 2pt + FH stability: `python stability_plot.py --tmin_max_stability --tpt_fh`
- Fig. 14+16 2pt + 3pt stability: `python stability_plot.py --tmin_max_stability --tpt_3pt`
- Fig. 15 [in the works]:
- Fig. 17+18 combined 2pt, 3pt, fh stability: `python stability_plot.py --tmin_max_stability --tpt_3pt_fh`
- Fig. 19 prior width stability: `python stability_plot.py --tmin_max_stability --prior_width`

Figures will be saved in a `new_plots` folder.

## Do the fit with database

```
fit_with_database.py

you can change parameters to do any fit you want, and choose to store results in the database or not.
```

## About code in module folder

```
fit.py : all fit functions and process about fits

p0.py : best_p0 to speed up fit process

plot.py : all plot functions

prepare_data.py : deal with the h5 data file, do average and output a dict

prior_setting.py : all prior models

```
