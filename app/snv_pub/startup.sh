#!/bin/bash
python manage.py createcachetable
python manage.py collectstatic --no-input
python manage.py makemigrations
python manage.py migrate
# gunicorn snv_test.wsgi --bind 0.0.0.0:8000 --timeout 600
python manage.py runserver 0.0.0.0:8000