"""Admin pages for ff_fit models

On default generates list view admins for all models
"""
from espressodb.base.admin import register_admins

register_admins("fh_db.ff_fit")
