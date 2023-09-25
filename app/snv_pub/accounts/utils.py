from allauth.account import app_settings
from allauth.account.models import EmailAddress, EmailConfirmation
from allauth.account.utils import get_adapter, user_email, user_pk_to_url_str
from datetime import timedelta
from django.contrib import messages
from django.utils.timezone import now
from logging import getLogger


logger_default = getLogger(__name__)


def send_email_confirmation(request, user, signup=False, email=None, logger=logger_default):
    """
    This is an improvement of allauth.account.utils.send_email_confirmation,
    which now accepts email as an argument.

    E-mail verification mails are sent:
    a) Explicitly: when a user signs up
    b) Implicitly: when a user attempts to log in using an unverified
    e-mail while EMAIL_VERIFICATION is mandatory.

    Especially in case of b), we want to limit the number of mails
    sent (consider a user retrying a few times), which is why there is
    a cooldown period before sending a new mail. This cooldown period
    can be configured in ACCOUNT_EMAIL_CONFIRMATION_COOLDOWN setting.
    """

    cooldown_period = timedelta(
        seconds=app_settings.EMAIL_CONFIRMATION_COOLDOWN
    )

    if not email:
        email = user_email(user)
    if email:
        try:
            email_address = EmailAddress.objects.get_for_user(user, email)
            if not email_address.verified:
                if app_settings.EMAIL_CONFIRMATION_HMAC:
                    send_email = True
                else:
                    send_email = not EmailConfirmation.objects.filter(
                        sent__gt=now() - cooldown_period,
                        email_address=email_address).exists()
                if send_email:
                    email_address.send_confirmation(request,
                                                    signup=signup)
            else:
                send_email = False
        except EmailAddress.DoesNotExist as e:
            logger.warning(e)
            pass

        if send_email:
            get_adapter(request).add_message(
                request,
                messages.INFO,
                'account/messages/'
                'email_confirmation_sent.txt',
                {'email': email})
    if signup:
        get_adapter(request).stash_user(request, user_pk_to_url_str(user))
