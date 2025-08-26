import os
import random
import logging
import base64
from io import BytesIO
from math import radians, sin, cos, sqrt, atan2
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import gpxpy
from django.shortcuts import render
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib import font_manager

logger = logging.getLogger(__name__)
matplotlib.use('Agg')

def distance(lat1, lon1, lat2, lon2):
    R = 6371.0
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return R * c

def reduce_data(moving_data, elevation_data, max_points=1000):
    if len(moving_data) > max_points:
        factor = len(moving_data) // max_points
        return moving_data[::factor], elevation_data[::factor]
    return moving_data, elevation_data

def reduce_gpx_data_based_on_length(gpx):
    points = gpx.tracks[0].segments[0].points
    total_points = len(points)
    step = 3 if total_points > 2000 else 2 if total_points > 1000 else 1
    reduced_points = points[::step]

    reduced_gpx = gpxpy.gpx.GPX()
    track = gpxpy.gpx.GPXTrack()
    segment = gpxpy.gpx.GPXTrackSegment(reduced_points)
    track.segments.append(segment)
    reduced_gpx.tracks.append(track)
    logger.info(f"Reduced points count: {len(reduced_points)}")
    return reduced_gpx

def meter_formatter(x, pos):
    return f'{int(x)}m'

def plot_elevation_profile(request):
    plt.clf()
    plt.cla()
    plt.close()

    if request.method == "POST":
        # Load and register HelveticaNeueLTPro font dynamically for this plot
        font_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts')

        # Register the font file HelveticaNeueLTPro.otf if it exists in fonts directory
        font_path = os.path.join(font_dir, 'HelveticaNeueLTPro.otf')
        if os.path.isfile(font_path):
            font_manager.fontManager.addfont(font_path)

        # Set matplotlib to use Helvetica Neue LT Pro as the font family
        plt.rcParams['font.family'] = 'Helvetica'

        gpx_file = request.FILES.get('gpx_file')
        aspect_ratio = float(request.POST.get('aspect', 2.5))
        line_color = request.POST.get('line_color', 'Blue')
        fill_color = request.POST.get('fill_color', '#E0EDD4')
        elements_color = request.POST.get('elements_color', '#52AE33')
        smoothness = int(request.POST.get('smoothness', 5))
        y_min = int(request.POST.get('y_min', 500))
        y_max = int(request.POST.get('y_max', 2000))
        grid_style = request.POST.get('grid_style', 'x')
        line_style = request.POST.get('line_style', 'solid')
        line_width = int(request.POST.get('line_width', 3))
        hide_box = request.POST.get('hide_box', False) == 'on'
        hide_grid = request.POST.get('hide_grid', False) == 'on'
        hide_labels = request.POST.get('hide_labels', False) == 'on'
        font_color = request.POST.get('font_color', '#666666')
        font_size = int(request.POST.get('font_size', 12))


        color_map = {
            "Blue": "#004a80",
            "Red": "#b91e2b",
            "Black": "#000000"
        }
        line_color = color_map.get(line_color, line_color)

        fill_color_rgb = mcolors.colorConverter.to_rgb(fill_color)
        color_rgb = mcolors.colorConverter.to_rgb(elements_color)

        hide = hide_box and hide_grid and hide_labels

        try:
            # --- Start: GPX-Datei bereinigen und parsen ---
            gpx_file.seek(0)
            gpx_content = gpx_file.read().decode('utf-8')
            gpx_content = gpx_content.replace('<![CDATA[', '')
            gpx_content = gpx_content.replace(']]>', '')
            gpx = gpxpy.parse(gpx_content)
            # --- Ende: GPX-Datei bereinigen und parsen ---

            # Verwende nur die Trackpoints für die Berechnung, da diese Höhenangaben haben
            # und den eigentlichen Trackverlauf darstellen.
            if not gpx.tracks or not gpx.tracks[0].segments:
                return render(request, 'elevation_app/index.html', {'error': 'Die GPX-Datei enthält keine auswertbaren Tracks.'})

            # Holen aller Punkte aus dem ersten Segment des ersten Tracks
            all_track_points = gpx.tracks[0].segments[0].points
            
            # Überprüfen, ob die Punkte eine Höhe haben
            if not any(p.elevation is not None for p in all_track_points):
                return render(request, 'elevation_app/index.html', {'error': 'Der Track enthält keine Höhenangaben.'})

            # Logging, um die Anzahl der Punkte zu überprüfen
            logger.info(f"Total track points: {len(all_track_points)}")
            
        except gpxpy.gpx.GPXXMLSyntaxException as e:
            logger.error(f"GPX XML Syntax Error: {e}")
            return render(request, 'elevation_app/index.html', {'error': 'Fehler: Die GPX-Datei hat ein ungültiges Format.'})
        except IndexError as e:
            logger.error(f"Index Error, likely no tracks in GPX: {e}")
            return render(request, 'elevation_app/index.html', {'error': 'Fehler: Die GPX-Datei enthält keine auswertbaren Tracks.'})
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            return render(request, 'elevation_app/index.html', {'error': 'Ein unerwarteter Fehler ist aufgetreten.'})

        # Jetzt Elevation und Moving Data berechnen
        elevation_data = [p.elevation for p in all_track_points]
        
        # Berechnung der kumulativen Distanz unter Verwendung der 3D-Distanz von gpxpy
        moving_data = [0]
        for i in range(1, len(all_track_points)):
            p1 = all_track_points[i-1]
            p2 = all_track_points[i]
            
            dist_3d = p1.distance_3d(p2)
            
            if dist_3d is not None:
                moving_data.append(moving_data[-1] + dist_3d / 1000) # In km
            else:
                # Dieser Fall sollte nun nicht mehr auftreten, aber als Fallback
                dist_2d = p1.distance_2d(p2)
                moving_data.append(moving_data[-1] + dist_2d / 1000)

        name = os.path.splitext(os.path.basename(gpx_file.name))[0]
        
        # Reduktion der Daten, falls nötig, um die Performance zu verbessern
        if len(moving_data) > 1000:
            moving_data, elevation_data = reduce_data(moving_data, elevation_data, max_points=1000)


        x = np.array(moving_data)
        y = np.array(elevation_data)

        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
        ax = plt.gca()

        if hide:
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_xlabel("Distance (km)", fontsize=12)
            ax.set_ylabel("", fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=10)

        if hide_grid:
            plt.grid(False)
        else:
            if grid_style == "both":
                plt.grid(True)
            elif grid_style == "x":
                plt.grid(True)
                ax.tick_params(axis='x', grid_linewidth=0)
            elif grid_style == "y":
                plt.grid(True)
                ax.tick_params(axis='y', grid_linewidth=0)
            else:
                plt.grid(False)

        plt.box(not hide_box)

        max_ticks = 9
        max_distance = moving_data[-1]

        if max_distance == 0:
            dash_line_moving = np.array([0])
        else:
            full_km_max = int(np.floor(max_distance))

            # Find the largest possible step that gives ≤ max_ticks (excluding 0, which we add manually)
            for step in range(1, full_km_max + 1):
                num_ticks = (full_km_max // step) + 1  # +1 for the 0 km tick
                if num_ticks <= max_ticks:
                    break

            # Create ticks from step to full_km_max
            dash_line_moving = np.arange(step, full_km_max + 1, step)

            # Prepend 0 km
            dash_line_moving = np.insert(dash_line_moving, 0, 0)

            # Append actual max_distance if not a full km
            if not np.isclose(max_distance, full_km_max):
                dash_line_moving = np.append(dash_line_moving, round(max_distance, 1))

        dash_line_elevation = [
            elevation_data[np.argmin(np.abs(moving_data - d))]
            for d in dash_line_moving
        ]



        y_max_plot = max(elevation_data) + 150 if max(elevation_data) > 1200 else 1200
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(xmin, xmax)
        ax.tick_params(axis='both', labelcolor="black", length=0, width=0)
        ax.grid(color="#9d9d9c")
        ax.set_box_aspect(1 / aspect_ratio)

        # x-axis integer ticks without commas
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        # y-axis with 'm' suffix
        ax.yaxis.set_major_formatter(FuncFormatter(meter_formatter))

        # Set y-axis tick label color (80% black) and font
        for label in ax.get_yticklabels():
            label.set_color("#666666")  # 80% black
            label.set_fontname("Helvetica")

        dash_line_moving_xticks = dash_line_moving
        ax.tick_params(axis='x', colors=elements_color, labelsize=10)
        plt.xticks(dash_line_moving_xticks)

        moving_data_array = np.array(moving_data)
        elevation_data_array = np.array(elevation_data)
        counts = {x: 0 for x in moving_data_array}
        for x in moving_data_array:
            counts[x] += 1
        epsilon = 1e-6
        new_moving_data = [
            x + random.uniform(-epsilon, epsilon) if counts[x] > 1 else x
            for x in moving_data_array
        ]
        new_moving_data_array = np.array(new_moving_data)

        new_moving_data_smooth = savgol_filter(new_moving_data_array, smoothness, 2)
        elevation_data_smooth = savgol_filter(elevation_data_array, smoothness, 2)

        f = interp1d(new_moving_data_array, elevation_data_smooth, kind='slinear', fill_value='extrapolate')
        xnew = np.linspace(new_moving_data_array.min(), new_moving_data_array.max(), len(new_moving_data) * 10)
        xnew = np.clip(xnew, moving_data_array.min(), moving_data_array.max())
        last_tick = dash_line_moving[-1]
        if last_tick > xnew.max():
            xnew = np.append(xnew, last_tick)

        # print("xnew min/max:", xnew.min(), xnew.max())
        # print("last_tick:", last_tick)
        ynew = f(xnew)

        def add_rectangles():
            rect_height = 15
            for i in range(len(dash_line_moving) - 1):
                color_fill = "#EFEFEF" if i % 2 == 0 else elements_color
                width = dash_line_moving[i + 1] - dash_line_moving[i]
                ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), width, rect_height, facecolor=color_fill, alpha=1, zorder=3))
            last_color = elements_color if len(dash_line_moving) % 2 == 0 else "#EFEFEF"
            last_width = max_distance - dash_line_moving[-1]
            if last_width > 0:
                ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), last_width, rect_height, facecolor=last_color, alpha=1, zorder=3))

        add_rectangles()

        plt.rcParams['svg.fonttype'] = 'none'
        plt.vlines(dash_line_moving, ymin=0, ymax=dash_line_elevation, colors=elements_color, linestyle='dashed', linewidth=2, zorder=2)
        plt.fill_between(xnew, ynew, color=fill_color, alpha=0.9, zorder=0.5)
        plt.plot(xnew, ynew, color=line_color, linestyle=line_style, linewidth=line_width, zorder=4)
        ax.set_title(name)

        buffer = BytesIO()
        plt.savefig(buffer, dpi=500, transparent=False, format='png')
        buffer.seek(0)
        graph = base64.b64encode(buffer.getvalue()).decode('utf-8')
        buffer.close()

        media_dir = 'media'
        os.makedirs(media_dir, exist_ok=True)

        svg_file_name = f"{name}_elevation_profile.svg"
        pdf_file_name = f"{name}_elevation_profile.pdf"

        svg_buffer = BytesIO()
        plt.savefig(svg_buffer, dpi=1000, transparent=False, format='svg', bbox_inches='tight', pad_inches=0.5)
        svg_buffer.seek(0)
        with open(os.path.join(media_dir, svg_file_name), 'wb') as f_svg:
            f_svg.write(svg_buffer.getvalue())
        svg_buffer.close()

        pdf_buffer = BytesIO()
        plt.savefig(pdf_buffer, dpi=1000, transparent=False, format='pdf', bbox_inches='tight', pad_inches=0.5)
        pdf_buffer.seek(0)
        with open(os.path.join(media_dir, pdf_file_name), 'wb') as f_pdf:
            f_pdf.write(pdf_buffer.getvalue())
        pdf_buffer.close()

        return render(request, 'elevation_app/index.html', {'graph': graph, 'svg_file_name': svg_file_name, 'pdf_file_name': pdf_file_name})

    return render(request, 'elevation_app/index.html')
