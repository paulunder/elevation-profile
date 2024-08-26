import random
from django.shortcuts import render
from io import BytesIO
import base64
from matplotlib import patches
import numpy as np
import gpxpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.signal import savgol_filter
import os
from math import radians, cos, sin, sqrt, atan2
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
import logging

logger = logging.getLogger(__name__)
def distance(lat1, lon1, lat2, lon2):
    # Radius of the Earth in kilometers
    R = 6371.0

    # Convert latitude and longitude from degrees to radians
    lat1 = radians(lat1)
    lon1 = radians(lon1)
    lat2 = radians(lat2)
    lon2 = radians(lon2)

    # Difference in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    # Distance in kilometers
    distance_km = R * c

    return distance_km

def reduce_data(moving_data, elevation_data, max_points=1000):
    if len(moving_data) > max_points:
        factor = len(moving_data) // max_points
        reduced_moving_data = moving_data[::factor]
        reduced_elevation_data = elevation_data[::factor]
        return reduced_moving_data, reduced_elevation_data
    else:
        return moving_data, elevation_data

def plot_elevation_profile(request):


    plt.clf()  # Clear the current plot
    plt.cla()  # Clear the current axes
    plt.close()

    if request.method == "POST":
        gpx_file = request.FILES.get('gpx_file')
        aspect_ratio = float(request.POST.get('aspect', 4))
        line_color = request.POST.get('line_color', 'Blue')
        fill_color = request.POST.get('fill_color', '#d4edfc')
        elements_color = request.POST.get('elements_color', '#004a80')
        smoothness = int(request.POST.get('smoothness', 5))
        y_min = int(request.POST.get('y_min', 500))
        y_max = int(request.POST.get('y_max', 1200))
        grid_style = request.POST.get('grid_style', 'x')
        line_style = request.POST.get('line_style', 'solid')
        line_width = int(request.POST.get('line_width', 3))
        hide_box = request.POST.get('hide_box', False) == 'on'
        hide_grid = request.POST.get('hide_grid', False) == 'on'
        hide_labels = request.POST.get('hide_labels', False) == 'on'
        color = elements_color


        if line_color == "Blue":
            line_color = "#004a80"
        elif line_color == "Red":
            line_color = "#b91e2b"
        elif line_color == "Black":
            line_color = "#000000"

        #y_label = values["-YLABEL-"]
        fill_color_rgb = mcolors.colorConverter.to_rgb(fill_color)
        color_rgb = mcolors.colorConverter.to_rgb(color)

        #setup the font
        #font = "Helvetica"


        #setup specific plot parameters
        fill_level = 50
        top_level = 70
        bottom_level = 50



        #make a boolean of values["-HIDE-"] to check if the box, ticks, title and labels should be hidden
        if hide_box and hide_grid and hide_labels:
            hide = True
        else:
            hide = False

        if hide_grid == True:
            hide_grid = True
        else:
            hide_grid = False

        if hide_box == True:
            hide_box = True
        else:
            hide_box = False



        # Load the GPX data
        gpx = gpxpy.parse(gpx_file)

        # Extract the elevation data
        elevation_data = [p.elevation for p in gpx.tracks[0].segments[0].points]

        # extract lat and lon data
        lat_data = [p.latitude for p in gpx.tracks[0].segments[0].points]
        lon_data = [p.longitude for p in gpx.tracks[0].segments[0].points]


        #define function named distance to calculate the distance between two points
        def distance(lat1, lon1, lat2, lon2):
            # approximate radius of earth in km
            R = 6373.0

            lat1 = np.radians(lat1)
            lon1 = np.radians(lon1)
            lat2 = np.radians(lat2)
            lon2 = np.radians(lon2)

            dlon = lon2 - lon1
            dlat = lat2 - lat1

            a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
            c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

            distance = R * c
            return distance

        #calculate the distance between the points
        moving_data = [0]
        for i in range(1, len(lat_data)):
            moving_data.append(moving_data[i-1] + distance(lat_data[i-1], lon_data[i-1], lat_data[i], lon_data[i]))

        #convert the moving_data from km to m
        #moving_data = [i * 1000 for i in moving_data]


        #extract the the file name of the input path
        name = os.path.basename(gpx_file.name)

        #delete the file extension from the file name
        name = os.path.splitext(name)[0]

        #import the aspect ratio and convert it to a float
        aspect = float(aspect_ratio)



        #a function loop to reduce the moving_data and the elevation_data by slicing it /2 and repeat it until it has at maximum 500 points
        def reduce_data(moving_data, elevation_data):
            #make a while loop to repeat the function until the length of the moving_data is less than 500
            while len(moving_data) > 800:
                #slice the moving_data and elevation_data by /2
                moving_data = moving_data[::2]
                elevation_data = elevation_data[::2]
            #return the reduced moving_data and elevation_data
            return moving_data, elevation_data


        #execute the function with an if statement to check if the length of the moving_data and the elevation_data is more than 500 points
        if len(moving_data) > 500:
            moving_data, elevation_data = reduce_data(moving_data, elevation_data)
        else:
            pass


            x = np.array(moving_data)
            y = np.array(elevation_data)


            xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
            ax = plt.gca()


            if hide == True:
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                ax.set_xlabel("Distance (km)")
                ax.set_ylabel("")
                ax.tick_params(axis='both', which='major', labelsize=10)

                #make a statement for -GRIDSTYLE- 
            if grid_style == "both":
                plt.grid(True)
            elif grid_style == "x":
                plt.grid(True)
                ax.tick_params(axis='x', grid_linewidth=0)
            elif grid_style == "y":
                plt.grid(True)
                ax.tick_params(axis='y', grid_linewidth=0)
            elif grid_style == "none":
                plt.grid(False)


            if hide_box == False:
                plt.box(True)
                # for tick in ax.get_xticklabels():
                #     tick.set_fontname(font)
                # for tick in ax.get_yticklabels():
                #     tick.set_fontname(font)
            else:
                plt.box(False)
                #ax.spines['left'].set_position('center')
                
                #ax.spines['left'].set_zorder(1000)
                #customise the ticks box itself
                # for tick in ax.get_xticklabels():
                #     tick.set_fontname(font)
                #     # tick.set_horizontalalignment("center")

                # for tick in ax.get_yticklabels():
                #     tick.set_fontname(font)
                #     # tick.set_horizontalalignment("center")
                    
            

            #make a new array and with these values 0.5km, 1km, 1.5km, 2km, 2.5km, 3km, 3.5km, 4km and so on until the maximum distance of the moving_data, make a statement if moving_data is above 5km then just add 1km, 2km, 3km, 4km and so on to the array
            if moving_data[-1] < 1:
                dash_line_moving = np.arange(0.1, moving_data[-1], 0.1)
            elif moving_data[-1] > 5:
                dash_line_moving = np.arange(1, moving_data[-1], 1)
            else:
                dash_line_moving = np.arange(0.5, moving_data[-1], 0.5)

            #make a statement - if the last value of moving_data is near (from 0 to 0.2) to the last value of dash_line_moving, then delete the last value of dash_line_moving
            if moving_data[-1] > 1:
                if moving_data[-1] - dash_line_moving[-1] < 0.2:
                    #create a rectangle of the last value of dash_line_moving
                    #ax.add_patch(patches.Rectangle((dash_line_moving[-1], ymin), 1, 15, facecolor=color, alpha=1, zorder=3))
                    dash_line_moving = np.delete(dash_line_moving, -1)
                    print("delete last value of dash_line_moving")
            else:
                pass

            dash_line_elevation = []
            for d in dash_line_moving:
            # find the index of the closest value in moving_data to d
                idx = np.argmin(np.abs(moving_data - d))
            # get the elevation value at that index from elevation_data
                e = elevation_data[idx]
            # append the elevation value to the dash_line_elevation list
                dash_line_elevation.append(e)

            if max(elevation_data) > 1200:
                y_max = max(elevation_data) + 150
            else:
                y_max = 1200


            ax.set_ylim(y_min, y_max)
            ax.set_xlim(xmin, xmax)
            print("moving_data: ")
            print(moving_data)
            print("---------------------")
            print("elevation_data: ")
            print(elevation_data)
            print("---------------------")
            print("Length of moving_data: ")
            print(len(moving_data))
            print("---------------------")
            print("Length of elevation_data:")
            print(len(elevation_data))

            #exclude the "-" from the tick labels and customise the grid lines
            ax.tick_params(axis='both', labelcolor="black", length=0, width=0)

            # ax.set_axisbelow(True)

            #set the gridline color
            ax.grid(color="#9d9d9c")

            ax.set_box_aspect(1/aspect)



            #FOR PLOTS WITH SMALL DISTANCE
            #add a last value to the dash_line_moving array, take the last value of the array and add + 0.5
            dash_line_moving_xticks = np.append(dash_line_moving, moving_data[-1])
            
                



            #make the xticks in clor = color
            ax.tick_params(axis='x', colors=color, labelsize=10)
            
            print("----------------")
            print("xtick locations: ")
            print(dash_line_moving_xticks)
            print("-----------------")
            print("dash line locations: ")
            print(dash_line_moving)
            print("-----------------")



            #smooth out the line with interp1d
            # Your original data
            moving_data_array = np.array(moving_data)
            moving_data_min = moving_data_array.min()
            moving_data_max = moving_data_array.max()
            elevation_data_array = np.array(elevation_data)

            # Create a dictionary to store the counts of each value
            counts = {}
            for x in moving_data_array:
                counts[x] = counts.get(x, 0) + 1

            # Add a small random noise to the duplicates
            epsilon = 1e-6 # You can adjust this value as needed
            new_moving_data = []
            for x in moving_data_array:
                if counts[x] > 1: # If x is a duplicate
                    x += random.uniform(-epsilon, epsilon) # Add some noise
                new_moving_data.append(x)

            # Convert the list back to an array
            new_moving_data_array = np.array(new_moving_data)
            new_moving_data_array_min = new_moving_data_array.min()
            new_moving_data_array_max = new_moving_data_array.max()

            #smooth out the new_moving_data_array and elevation_data_array using savgol_filter
            new_moving_data_array_smooth = savgol_filter(new_moving_data_array, smoothness, 2)
            elevation_data_array_smooth = savgol_filter(elevation_data_array, smoothness, 2)

            # Smooth out the line using interp1d
            f = interp1d(new_moving_data_array, elevation_data_array_smooth, kind='slinear')
            xnew = np.linspace(new_moving_data_array_min, new_moving_data_array_max, len(new_moving_data)*10)
            # Clip the values of xnew to the range of moving_data
            xnew = np.clip(xnew, moving_data_min, moving_data_max)
            ynew = f(xnew)




            plt.xticks(dash_line_moving_xticks)

            if moving_data[-1] > 5:
                print("create rectangles - track is above 5km")
                #create a rectangle at y y_min to 460 and x 0 to 0.5
                ax.add_patch(patches.Rectangle((0, y_min), 1, 15, facecolor=color, alpha=1, zorder=3))

                for i in range(0, len(dash_line_moving)-1):
                    if i % 2 == 0:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
                    else:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=color, alpha=1, zorder=3))

                #create the last rectangle 
                if int(len(dash_line_moving)) % 2 == 0:
                    ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor=color, alpha=1, zorder=3))
                else:
                    ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))

            elif moving_data[-1] < 1:
                print("create rectangles - track is under 1km")
                #create a rectangle at y y_min to 460 and x from 0 to first value in dash_line_moving
                ax.add_patch(patches.Rectangle((0, y_min), dash_line_moving[0], 15, facecolor=color, alpha=1, zorder=3))

                #create a loop that will repeat the above code for every value in dash_line_moving, stop the loop at the last value in dash_line_moving, 
                #also create alway the first rectangle from 0 to the first value in dash_line_moving

                for i in range(0, len(dash_line_moving)-1):
                    if i % 2 == 0:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
                    else:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=color, alpha=1, zorder=3))
                
                #create a rectangle at y y_min to 460 and x from prelast value in dash_line_moving to last value in dash_line_moving, choose the opposite color based on the last rectangle created
                if i % 2 == 0:
                    ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), 1, 15, facecolor=color, alpha=1, zorder=3))


            #make a statement elif, if moving_data[-1] is between 1 and 5
            elif moving_data[-1] > 1 and moving_data[-1] < 5:
                print("rectangles - track is between 1km and 5km")

                #create a loop that will repeat the above code for every value in dash_line_moving, stop the loop at the last value in dash_line_moving
                for i in range(0, len(dash_line_moving)-1):
                    if i % 2 == 0:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
                    else:
                        ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=color, alpha=1, zorder=3))

                #create a rectangle at y y_min to 460 and x 0 to 0.5
                ax.add_patch(patches.Rectangle((0, y_min), 0.5, 15, facecolor=color, alpha=1, zorder=3))

                #create a rectangle at y y_min to 460 and x from prelast value in dash_line_moving to last value in dash_line_moving, choose the opposite color based on the last rectangle created
                if int(len(dash_line_moving)) % 2 == 0:
                    ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor=color, alpha=1, zorder=3))
                else:
                    ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))

            #fix font bug in svg files
            plt.rcParams['svg.fonttype'] = 'none'

            #make a dashed line vertical which ends at the line (elevation data), the line should be plotted for dash_line_moving
            plt.vlines(dash_line_moving, 0, dash_line_elevation, color=color, linestyle='dashed',linewidth=2, zorder=2)
            
            #make a field with no gradient
            plt.fill_between(xnew, ynew, color=fill_color, zorder=0.5, alpha=0.9)

            #plot the data
            plt.plot(xnew, ynew, color=line_color, linestyle=line_style, linewidth=line_width, zorder=4)

            ax.set_title(name)

            print("----------------")
            print("y_min: ")
            print(y_min)
            print("----------------")
            print("y_max: ")
            print(y_max)
            print("----------------")


        # Save the plot to a BytesIO object as PNG
            buffer = BytesIO()
            plt.savefig(buffer, dpi=500, transparent=False, format='png')
            buffer.seek(0)
            image_png = buffer.getvalue()
            graph = base64.b64encode(image_png).decode('utf-8')
            buffer.close()

            # Save the plot as SVG
            svg_buffer = BytesIO()
            plt.savefig(svg_buffer, dpi=1000, transparent=False, format='svg', bbox_inches='tight', pad_inches=0.5)
            svg_buffer.seek(0)
            svg_file_name = name + '_elevation_profile.svg'

            #save the plot as pdf
            pdf_buffer = BytesIO()
            plt.savefig(pdf_buffer, dpi=1000, transparent=False, format='pdf', bbox_inches='tight', pad_inches=0.5)
            pdf_buffer.seek(0)
            pdf_file_name = name + '_elevation_profile.pdf'

            # Ensure the 'media' directory exists
            media_dir = 'media'
            if not os.path.exists(media_dir):
                os.makedirs(media_dir)

            with open(os.path.join(media_dir, svg_file_name), 'wb') as f:
                f.write(svg_buffer.getvalue())
            svg_buffer.close()

            with open(os.path.join(media_dir, pdf_file_name), 'wb') as f:
                f.write(pdf_buffer.getvalue())
            pdf_buffer.close()

            return render(request, 'elevation_app/index.html', {'graph': graph, 'svg_file_name': svg_file_name})

    return render(request, 'elevation_app/index.html')


    # fig1 = plt.gcf()


    #     #make fig bigger, save it and show it
    #     fig1.set_size_inches(16, 9)
    #     fig1.savefig(values["-OUT-"], dpi=1000, transparent=False, format='svg')

    #     print("File saved: ")
    #     print(values["-OUT-"])

    #     plt.show()


############################################################################################################
    #     # Load the GPX data
    #     gpx = gpxpy.parse(gpx_file)

    #     # Extract the elevation and distance data
    #     elevation_data = [p.elevation for p in gpx.tracks[0].segments[0].points]
    #     lat_data = [p.latitude for p in gpx.tracks[0].segments[0].points]
    #     lon_data = [p.longitude for p in gpx.tracks[0].segments[0].points]

    #     # Calculate the moving data (distances between points)
    #     moving_data = [0]
    #     for i in range(1, len(lat_data)):
    #         moving_data.append(moving_data[i-1] + distance(lat_data[i-1], lon_data[i-1], lat_data[i], lon_data[i]))

    #     # Reduce the number of data points if necessary
    #     moving_data, elevation_data = reduce_data(moving_data, elevation_data)

    #     # Apply smoothing
    #     moving_data_array = np.array(moving_data)
    #     elevation_data_array = np.array(elevation_data)
    #     new_moving_data_array_smooth = savgol_filter(moving_data_array, smoothness, 2)
    #     elevation_data_array_smooth = savgol_filter(elevation_data_array, smoothness, 2)

    #     # Use a standard font family
    #     plt.rcParams['font.family'] = 'DejaVu Sans'
    #     plt.rcParams['svg.fonttype'] = 'none'  # Embed font in SVG

    #     y_min = min(elevation_data_array_smooth) - 10  # Optional buffer
    #     y_max = max(elevation_data_array_smooth) + 10  # Optional buffer

    #     # Create the plot
    #     fig, ax = plt.subplots()

    #     # Plot the smoothed data
    #     ax.plot(new_moving_data_array_smooth, elevation_data_array_smooth, color=line_color, linestyle=line_style, linewidth=line_width)

    #     # Set limits
    #     ax.set_xlim(min(new_moving_data_array_smooth), max(new_moving_data_array_smooth))
    #     ax.set_ylim(y_min, y_max)

    #     # Add patches based on the distance range
    #     if moving_data[-1] < 1:
    #         dash_line_moving = np.arange(0.1, moving_data[-1], 0.1)
    #     elif moving_data[-1] > 5:
    #         dash_line_moving = np.arange(1, moving_data[-1], 1)
    #     else:
    #         dash_line_moving = np.arange(0.5, moving_data[-1], 0.5)

    #     if moving_data[-1] > 1:
    #         if moving_data[-1] - dash_line_moving[-1] < 0.2:
    #             dash_line_moving = np.delete(dash_line_moving, -1)

    #     dash_line_elevation = []
    #     for d in dash_line_moving:
    #         idx = np.argmin(np.abs(moving_data - d))
    #         e = elevation_data[idx]
    #         dash_line_elevation.append(e)

    #     if moving_data[-1] > 5:
    #         ax.add_patch(patches.Rectangle((0, y_min), 1, 15, facecolor=elements_color, alpha=1, zorder=3))

    #         for i in range(0, len(dash_line_moving)-1):
    #             if i % 2 == 0:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
    #             else:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=elements_color, alpha=1, zorder=3))

    #         if int(len(dash_line_moving)) % 2 == 0:
    #             ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor=elements_color, alpha=1, zorder=3))
    #         else:
    #             ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))

    #     elif moving_data[-1] < 1:
    #         ax.add_patch(patches.Rectangle((0, y_min), dash_line_moving[0], 15, facecolor=fill_color, alpha=1, zorder=3))

    #         for i in range(0, len(dash_line_moving)-1):
    #             if i % 2 == 0:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
    #             else:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=elements_color, alpha=1, zorder=3))

    #         if i % 2 == 0:
    #             ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), 1, 15, facecolor=elements_color, alpha=1, zorder=3))

    #     elif moving_data[-1] > 1 and moving_data[-1] < 5:
    #         for i in range(0, len(dash_line_moving)-1):
    #             if i % 2 == 0:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))
    #             else:
    #                 ax.add_patch(patches.Rectangle((dash_line_moving[i], y_min), dash_line_moving[i+1], 15, facecolor=elements_color, alpha=1, zorder=3))

    #         ax.add_patch(patches.Rectangle((0, y_min), 0.5, 15, facecolor=elements_color, alpha=1, zorder=3))

    #         if int(len(dash_line_moving)) % 2 == 0:
    #             ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor=elements_color, alpha=1, zorder=3))
    #         else:
    #             ax.add_patch(patches.Rectangle((dash_line_moving[-1], y_min), moving_data[-1], 15, facecolor="#EFEFEF", alpha=1, zorder=3))

    #     # Set labels, grid, and box aspect
    #     if not hide_labels:
    #         ax.set_xlabel("Distance (km)")
    #         ax.set_ylabel("Elevation (m)")
    #         ax.set_title("Elevation Profile")

    #     if hide_box:
    #         ax.spines['top'].set_visible(False)
    #         ax.spines['right'].set_visible(False)
    #         ax.spines['left'].set_visible(False)
    #         ax.spines['bottom'].set_visible(False)
    #         ax.get_xaxis().set_ticks([])
    #         ax.get_yaxis().set_ticks([])

    #     if hide_grid:
    #         ax.grid(False)
    #     else:
    #         ax.grid(which=grid_style)

    #     ax.set_box_aspect(1/aspect_ratio)

    #     # Save the plot to a BytesIO object as PNG
    #     buffer = BytesIO()
    #     plt.savefig(buffer, format='png')
    #     buffer.seek(0)
    #     image_png = buffer.getvalue()
    #     graph = base64.b64encode(image_png).decode('utf-8')
    #     buffer.close()

    #     # Save the plot as SVG
    #     svg_buffer = BytesIO()
    #     plt.savefig(svg_buffer, format='svg')
    #     svg_buffer.seek(0)
    #     svg_file_name = 'elevation_profile.svg'

    #     # Ensure the 'media' directory exists
    #     media_dir = 'media'
    #     if not os.path.exists(media_dir):
    #         os.makedirs(media_dir)

    #     with open(os.path.join(media_dir, svg_file_name), 'wb') as f:
    #         f.write(svg_buffer.getvalue())
    #     svg_buffer.close()

    #     return render(request, 'elevation_app/index.html', {'graph': graph, 'svg_file_name': svg_file_name})

    # return render(request, 'elevation_app/index.html')
