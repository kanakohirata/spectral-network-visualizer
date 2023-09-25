from allauth.account import signals
from allauth.account.forms import AddEmailForm
from allauth.account.models import EmailAddress
from allauth.account.utils import filter_users_by_email, get_adapter
from allauth.account.views import AjaxCapableProcessFormViewMixin, EmailView, _ajax_response
from allauth.decorators import rate_limit
from allauth.utils import get_form_class, get_request_param
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.sites.shortcuts import get_current_site
from django.core.mail import send_mail
from django.dispatch import receiver
from django.http import HttpResponseRedirect
from django.shortcuts import render, redirect
from django.template.loader import render_to_string
from django.urls import reverse_lazy
from django.utils.decorators import method_decorator
from django.views import generic
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from logging import getLogger

from accounts.forms.user import DeleteAccountFrom
from accounts.models import CustomUser
from accounts.utils import send_email_confirmation

logger = getLogger(__name__)


class UserProfileView(generic.TemplateView):
    http_method_names = ['get']
    template_name = 'setting/profile.html'
    success_url = reverse_lazy('accounts:profile')
    email_add_form_class = AddEmailForm
    delete_account_form_class = DeleteAccountFrom

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        primary_email = EmailAddress.objects.get(user=self.request.user, primary=True)
        other_emails = EmailAddress.objects.filter(user=self.request.user, primary=False)
        logger.debug(other_emails)
        context['primary_email'] = primary_email
        context['other_emails'] = other_emails
        context['email_add_form'] = self.email_add_form_class()
        context['delete_account_form'] = self.delete_account_form_class()

        return context


user_profile = login_required(UserProfileView.as_view())


class CustomEmailView(EmailView):
    http_method_names = ['post', 'get']
    template_name = 'setting/profile.html'
    success_url = reverse_lazy('accounts:profile')
    email_add_form_class = AddEmailForm
    delete_account_form_class = DeleteAccountFrom

    def get(self):
        redirect('accounts:profile')

    def form_invalid(self, form):
        """If the form is invalid, render the invalid form."""
        return self.render_to_response(self.get_context_data(email_add_form=form))

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        primary_email = EmailAddress.objects.get(user=self.request.user, primary=True)
        other_emails = EmailAddress.objects.filter(user=self.request.user, primary=False)
        logger.debug(other_emails)
        context['primary_email'] = primary_email
        context['other_emails'] = other_emails
        context['delete_account_form'] = self.delete_account_form_class()

        if kwargs.get('email_add_form'):
            logger.warning(kwargs)
        else:
            context['email_add_form'] = self.email_add_form_class()

        return context


@method_decorator(rate_limit(action="manage_email"), name="dispatch")
class SendEmailVerificationView(CustomEmailView):

    def post(self, request, *args, **kwargs):
        if request.POST.get("email"):
            res = self._action_send(request)
            res = res or HttpResponseRedirect(self.get_success_url())
            # Given that we bypassed AjaxCapableProcessFormViewMixin,
            # we'll have to call invoke it manually...
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        else:
            # No email address selected
            res = HttpResponseRedirect(self.success_url)
            res = _ajax_response(request, res, data=self._get_ajax_data_if())

        return res


email_verify = login_required(SendEmailVerificationView.as_view())


@method_decorator(rate_limit(action="manage_email"), name="dispatch")
class RemoveEmailView(CustomEmailView):

    def post(self, request):
        if request.POST.get("email"):
            res = self._action_remove(request)
            res = res or HttpResponseRedirect(self.get_success_url())
            # Given that we bypassed AjaxCapableProcessFormViewMixin,
            # we'll have to call invoke it manually...
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        else:
            # No email address selected
            res = HttpResponseRedirect(self.success_url)
            res = _ajax_response(request, res, data=self._get_ajax_data_if())

        return res


email_remove = login_required(RemoveEmailView.as_view())


@method_decorator(rate_limit(action="manage_email"), name="dispatch")
class ChangePrimaryEmailView(CustomEmailView):

    def post(self, request):
        if request.POST.get("email"):
            to_email_address = EmailAddress.objects.get(user=request.user,
                                                        email=request.POST.get("email"))

            logger.warning(to_email_address.email)
            if to_email_address.primary:
                logger.warning(f'{to_email_address.email} is already primary email.')
                return redirect('accounts:profile')

            res = self._action_primary(request)
            res = res or HttpResponseRedirect(self.get_success_url())
            # Given that we bypassed AjaxCapableProcessFormViewMixin,
            # we'll have to call invoke it manually...
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        else:
            # No email address selected
            res = HttpResponseRedirect(self.success_url)
            res = _ajax_response(request, res, data=self._get_ajax_data_if())

        return res


email_change_primary = login_required(ChangePrimaryEmailView.as_view())


@method_decorator(rate_limit(action="manage_email"), name="dispatch")
class AddEmailView(CustomEmailView):

    def form_invalid(self, form):
        """If the form is invalid, render the invalid form."""
        return self.render_to_response(self.get_context_data(email_add_form=form))

    def post(self, request, *args, **kwargs):
        logger.debug(f'@@@@@ 0: {request.POST}')
        if "action_add" in request.POST:
            res = super(EmailView, self).post(request, *args, **kwargs)
        else:
            # No email address selected
            res = HttpResponseRedirect(self.success_url)
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        return res


email_add = login_required(AddEmailView.as_view())


@method_decorator(csrf_exempt, name='dispatch')
class DeleteAccountView(generic.FormView):
    form_class = DeleteAccountFrom
    template_name = 'setting/profile.html'
    success_url = 'visualizer:spectral_network'
    http_method_names = ['post', 'get']
    email_add_form_class = AddEmailForm
    delete_account_form_class = DeleteAccountFrom

    def get(self, request):
        return redirect('accounts:profile')

    @method_decorator(csrf_protect)
    def post(self, request, *args, **kwargs):
        res = super().post(request, *args, **kwargs)
        return res

    def form_valid(self, form):
        messages.success(self.request, 'Delete your account.')
        self.request.user.delete()
        return redirect(self.get_success_url())

    def form_invalid(self, form):
        messages.error(self.request, 'Invalid form!')
        res = self.render_to_response(self.get_context_data(delete_account_form=form,
                                                            delete_account_form_opened='open'))

        return res

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        primary_email = EmailAddress.objects.get(user=self.request.user, primary=True)
        other_emails = EmailAddress.objects.filter(user=self.request.user, primary=False)
        logger.debug(other_emails)
        context['primary_email'] = primary_email
        context['other_emails'] = other_emails
        context['email_add_form'] = self.email_add_form_class()

        if kwargs.get('delete_account_form'):
            pass
        else:
            context['delete_account_form'] = self.delete_account_form_class()

        return context


delete_account = login_required(DeleteAccountView.as_view())

