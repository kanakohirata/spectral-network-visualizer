from allauth.account import views
from allauth.account.views import _ajax_response
from django.contrib.auth.decorators import login_required
from django.http import HttpResponseRedirect
from django.shortcuts import render, redirect
from logging import getLogger


logger = getLogger(__name__)


class CustomAccountInactiveView(views.AccountInactiveView):
    def get(self, request, *args, **kwargs):
        if request.user.is_active:
            return redirect('spectrum:index')

        return super().get(request, *args, **kwargs)


account_inactive = CustomAccountInactiveView.as_view()


class EmailView(views.EmailView):
    """
    The post() function does not work correctly because the request.POST does not contain
    action types such as "action_add", "action_send", "action_remove" and "action_primary".
    In order to make it work, a javascript code is needed to get these values .
    """
    def post(self, request, *args, **kwargs):
        logger.debug('@@@@@ 0')
        logger.debug(f'request.POST: {request.POST}')
        res = None
        if "action_add" in request.POST:
            logger.debug('@@@@@1')
            form = self.form_class(request.user, **self.get_form_kwargs())
            if form.is_valid():
                logger.debug('@@@@@2 form is valid')
                res = self.form_valid(form)
            else:
                logger.debug('@@@@@2 form is invalid')
                res = self.form_invalid(form)
        elif request.POST.get("email"):
            if "action_send" in request.POST:
                res = self._action_send(request)
            elif "action_remove" in request.POST:
                res = self._action_remove(request)
            elif "action_primary" in request.POST:
                res = self._action_primary(request)
            res = res or HttpResponseRedirect(self.success_url)
            # Given that we bypassed AjaxCapableProcessFormViewMixin,
            # we'll have to call invoke it manually...
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        else:
            # No email address selected
            res = HttpResponseRedirect(self.success_url)
            res = _ajax_response(request, res, data=self._get_ajax_data_if())
        return res


class CustomEmailView(views.EmailView):
    http_method_names = ['get']

    def get(self, request, *args, **kwargs):
        return redirect('accounts:profile')


email = login_required(CustomEmailView.as_view())
