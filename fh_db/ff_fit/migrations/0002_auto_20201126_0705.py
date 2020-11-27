# Generated by Django 3.0.6 on 2020-11-26 07:05

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ff_fit', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ff_fit',
            name='id_num',
            field=models.PositiveSmallIntegerField(help_text='used to divide and select', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt2_tmax',
            field=models.PositiveSmallIntegerField(help_text='Maximum of t in 2pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt2_tmin',
            field=models.PositiveSmallIntegerField(help_text='Minimum of t in 2pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt3_A3_tsep_max',
            field=models.PositiveSmallIntegerField(help_text='Maximum of tsep in 3pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt3_A3_tsep_min',
            field=models.PositiveSmallIntegerField(help_text='Minimum of tsep in 3pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt3_V4_tsep_max',
            field=models.PositiveSmallIntegerField(help_text='Maximum of tsep in 3pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='pt3_V4_tsep_min',
            field=models.PositiveSmallIntegerField(help_text='Minimum of tsep in 3pt fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='sum_A3_tsep_max',
            field=models.PositiveSmallIntegerField(help_text='Maximum of tsep in sum tau fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='sum_A3_tsep_min',
            field=models.PositiveSmallIntegerField(help_text='Minimum of tsep in sum tau fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='sum_V4_tsep_max',
            field=models.PositiveSmallIntegerField(help_text='Maximum of tsep in sum tau fit', null=True),
        ),
        migrations.AlterField(
            model_name='ff_fit',
            name='sum_V4_tsep_min',
            field=models.PositiveSmallIntegerField(help_text='Minimum of tsep in sum tau fit', null=True),
        ),
    ]