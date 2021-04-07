# READ ME

This project is used to fit a09 data and produce all plots in fh comparison paper.

## Prepare

In order to use database:

```
1. run "pip install espressodb"

2. change the "NAME" in "db-config.yaml" file to the correct path of "fh_db-db.sqlite" file

3. run "python manage.py makemigrations"

4. run "python manage.py migrate"
```

## Figures in the paper

```
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