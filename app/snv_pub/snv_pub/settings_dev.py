from .settings_common import *

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = False

ALLOWED_HOSTS = ['*']

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'whitenoise.runserver_nostatic',
    'django.contrib.staticfiles',
    'accounts.apps.AccountsConfig',
    'visualizer.apps.VisualizerConfig',

    # django-allauth requires the following 4 apps
    'django.contrib.sites',
    'allauth',
    'allauth.account',
    'allauth.socialaccount',

    # Use 'debug_toolbar' only when DEBUG = True
    'debug_toolbar',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'whitenoise.middleware.WhiteNoiseMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',

    # debug_toolbar needs the following.
    'debug_toolbar.middleware.DebugToolbarMiddleware',
]

# debug_toolbar needs the following.
INTERNAL_IPS = ['127.0.0.1']


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'OPTIONS': {
            'options': f'-c search_path={env("DB_SCHEMA")},public',
        },
        'NAME': 'snv_pub',
        'USER': env('DB_USER'),
        'PASSWORD': env('DB_PASSWORD'),
        'HOST': env('DB_HOST'),
        'PORT': env('DB_PORT'),
        'TEST': {
            'NAME': 'my_test',
        }
    }
}

if not os.path.isdir(f'{BASE_DIR}/visualizer/processing/multilayer_3d_ms_network/logs'):
    os.makedirs(f'{BASE_DIR}/visualizer/processing/multilayer_3d_ms_network/logs')

if not os.path.isdir(f'{BASE_DIR}/report'):
    os.makedirs(f'{BASE_DIR}/report')


LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,

    # Logger
    'loggers': {
        # Logger for Django
        'django': {
            'handlers': ['console'],
            'level': 'INFO',
        },
        # Logger for visualizer app
        'visualizer': {
            'handlers': ['console'],
            'level': 'DEBUG',
        },
        # Logger for multilayer_3d_ms_network package
        'visualizer.processing.multilayer_3d_ms_network': {
            'handlers': ['console', 'file_multilayer_3d_ms_network'],
            'level': 'WARNING',
            'propagate': False,
        },
        'visualizer.report': {
            'handlers': ['visualizer_report'],
            'level': 'INFO',
            'propagate': False,
        }
    },

    # Handler
    'handlers': {
        'console': {
            'level': 'DEBUG',
            'class': 'logging.StreamHandler',
            'formatter': 'dev'
        },
        'file_multilayer_3d_ms_network': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'formatter': 'dev',
            'filename': 'visualizer/processing/multilayer_3d_ms_network/logs/info.log',
        },
        'visualizer_report': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'formatter': 'report',
            'filename': 'report/visualizer.log',
        }
    },

    # Formatter
    'formatters': {
        'dev': {
            'format': '\t'.join([
                '%(asctime)s',
                '[%(levelname)s]',
                '%(name)s(Line:%(lineno)d)',
                '%(message)s'
            ])
        },
        'report': {
            'format': '%(message)s'
        },
    }
}

STATIC_ROOT = os.path.join(BASE_DIR, 'staticall')
if not os.path.isdir(STATIC_ROOT):
    os.makedirs(STATIC_ROOT)

FILE_UPLOAD_TEMP_DIR = os.path.join(BASE_DIR, 'tmp')
if not os.path.isdir(FILE_UPLOAD_TEMP_DIR):
    os.makedirs(FILE_UPLOAD_TEMP_DIR)
