# READ ME

This project is used to fit a09 data and produce all plots in fh comparison paper.

## Prepare

In order to use database:

```
1. run "pip install espressodb"

2. change the "NAME" in "db-config.yaml" file to the correct path of "fh_db-db1.sqlite" file (this is not a completed database)

3. run "python manage.py makemigrations"

4. run "python manage.py migrate"
```

## Figures in the paper

```
you can make these plots in the paper with corresponding code

Fig.1: fit_on_data_plot.py
Fig.3: stability_plot.py
Fig.4: excited_state.py
Fig.6: stability_plot.py
Fig.7: stability_plot.py
Fig.8: stability_plot.py
Fig.10: fit_on_data_plot.py
Fig.11: stability_plot.py
Fig.12: stability_plot.py
Fig.13: stability_plot.py
Fig.15: stability_plot.py
Fig.16: stability_plot.py
Fig.17: stability_plot.py
Fig.18: stability_plot.py

```

## Do the fit with database

```
fit_with_database.py

you can change parameters to do any fit you want, and choose to store results in the database or not.
```