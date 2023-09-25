from allauth.account import app_settings
from allauth.account.forms import UserForm
from allauth.account.models import EmailAddress
from allauth.account.utils import filter_users_by_email, get_adapter
from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext, gettext_lazy as _, pgettext
from logging import getLogger

from accounts.models import CustomUser

logger = getLogger(__name__)


class DeleteAccountFrom(UserForm):

    confirmation = forms.CharField(max_length=10, required=True)

    def clean(self):
        cleaned_data = super().clean()
        confirmation = cleaned_data.get('confirmation')

        if confirmation != 'delete':
            raise ValidationError(
                {'confirmation': 'Please type "delete" to delete your account.'}
            )
