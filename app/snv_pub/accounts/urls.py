from django.urls import path
from accounts import views

app_name = 'accounts'

urlpatterns = [
    path('profile/', views.user_profile, name='profile'),
    path('email/add', views.email_add, name='email_add'),
    path('email/verify', views.email_verify, name='email_verify'),
    path('email/remove', views.email_remove, name='email_remove'),
    path('email/change-primary', views.email_change_primary, name='email_change_primary'),
    path('delete-account', views.delete_account, name='delete_account'),
]
