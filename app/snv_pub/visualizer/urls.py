from django.urls import path
from . import views

app_name = 'visualizer'

urlpatterns = [
    path('', views.SpectralNetworkView.as_view(), name='spectral_network'),
    path('visualizer/network-data/', views.GetNetworkData.as_view(), name='get_network_data'),
]
