"""Models of ff_fit
"""

# Note: if you want your models to use espressodb features, they must inherit from Base

from django.db import models
from espressodb.base.models import Base
class FF_fit(Base):
    data_file_name = models.TextField(
        null=False,
        blank=False,
        help_text="Store the name of data file",
    )
    include_2pt = models.BooleanField(
        null=True, help_text="Whether include 2pt part in fitting"
    )
    include_3pt = models.BooleanField(
        null=False, help_text="Whether include 3pt part in fitting"
    )
    include_sum = models.BooleanField(
        null=False, help_text="Whether include sum part in fitting"
    )
    data_and_results = models.TextField(
        null=False, help_text="Use gv.dump to package data, prama, chi2 and dof, and turn into hexcode"
    )
    prior_hexcode = models.TextField(
        null=False, help_text="put all priors into a string with fixed order and turn it into hexcode"
    )
    pt2_nstates = models.PositiveSmallIntegerField(
        null=False, help_text="Number of states in 2pt fit"
    )
    pt2_tmin = models.PositiveSmallIntegerField(
        null=True, help_text="Minimum of t in 2pt fit"
    )
    pt2_tmax = models.PositiveSmallIntegerField(
        null=True, help_text="Maximum of t in 2pt fit"
    )
    pt3_nstates = models.PositiveSmallIntegerField(
        null=False, help_text="Number of states in 3pt fit"
    )
    pt3_A3_tsep_min = models.PositiveSmallIntegerField(
        null=True, help_text="Minimum of tsep in 3pt fit"
    )
    pt3_A3_tsep_max = models.PositiveSmallIntegerField(
        null=True, help_text="Maximum of tsep in 3pt fit"
    )
    pt3_V4_tsep_min = models.PositiveSmallIntegerField(
        null=True, help_text="Minimum of tsep in 3pt fit"
    )
    pt3_V4_tsep_max = models.PositiveSmallIntegerField(
        null=True, help_text="Maximum of tsep in 3pt fit"
    )
    pt3_tau_dict = models.TextField(
        null=True, help_text="Use gv.dump to package tau dict of 3pt, and turn into hexcode"
    )
    sum_nstates = models.PositiveSmallIntegerField(
        null=False, help_text="Number of states in sum tau fit"
    )
    sum_A3_tsep_min = models.PositiveSmallIntegerField(
        null=True, help_text="Minimum of tsep in sum tau fit"
    )
    sum_A3_tsep_max = models.PositiveSmallIntegerField(
        null=True, help_text="Maximum of tsep in sum tau fit"
    )
    sum_V4_tsep_min = models.PositiveSmallIntegerField(
        null=True, help_text="Minimum of tsep in sum tau fit"
    )
    sum_V4_tsep_max = models.PositiveSmallIntegerField(
        null=True, help_text="Maximum of tsep in sum tau fit"
    )
    sum_tau_cut = models.PositiveSmallIntegerField(
        null=False, help_text="Start of sum in sum tau fit"
    )
    E0 = models.FloatField(
        null=False, help_text="fit result of E0"
    )
    E0_err = models.FloatField(
        null=False, help_text="fit result of E0 err"
    )
    z0 = models.FloatField(
        null=False, help_text="fit result of z0"
    )
    z0_err = models.FloatField(
        null=False, help_text="fit result of z0 err"
    )
    z0_ps = models.FloatField(
        null=True, help_text="fit result of z0"
    )
    z0_ps_err = models.FloatField(
        null=True, help_text="fit result of z0 err"
    )
    A300 = models.FloatField(
        null=False, help_text="fit result of A300"
    )
    A300_err = models.FloatField(
        null=False, help_text="fit result of A300 err"
    )
    V400 = models.FloatField(
        null=False, help_text="fit result of V400"
    )
    V400_err = models.FloatField(
        null=False, help_text="fit result of V400 err"
    )
    Q_value = models.FloatField(
        null=False, help_text="Q value of fitting"
    )
    log_GBF = models.FloatField(
        null=False, help_text="logGBF of fitting"
    )
    id_num = models.PositiveSmallIntegerField(
        null=True, help_text="used to divide and select"
    )
    fit_type = models.TextField(
        null=False,
        help_text="separate or continue",
    )
    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["data_file_name", 'fit_type', 'include_2pt', 'include_3pt', 'include_sum', 'data_and_results', 'prior_hexcode','pt2_nstates','pt2_tmin','pt2_tmax',     'pt3_nstates','pt3_A3_tsep_min','pt3_A3_tsep_max','pt3_V4_tsep_min',                     'pt3_V4_tsep_max','sum_nstates',
                'pt3_tau_dict','sum_A3_tsep_min','sum_A3_tsep_max','sum_V4_tsep_min','sum_V4_tsep_max',                       'sum_tau_cut','id_num'], name="unique_fit"
            )
        ]