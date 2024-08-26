from django.urls import path
from .views import plot_elevation_profile

urlpatterns = [
    path('', plot_elevation_profile, name='plot_elevation_profile'),
]
