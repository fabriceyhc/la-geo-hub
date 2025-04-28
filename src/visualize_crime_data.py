# visualize_crime_data.py
import geopandas as gpd
import matplotlib.pyplot as plt
import contextily as ctx
import argparse
import os
import pandas as pd # Import pandas for numeric check

# --- Configuration ---
# Approximate bounding box for Los Angeles (min Lon, min Lat, max Lon, max Lat)
LA_BOUNDS = (-118.9, 33.7, -117.6, 34.35)
# Directory to save output images
IMG_OUTPUT_DIR = "imgs"
# Default field for color-coding
DEFAULT_COLOR_FIELD = "Crimes_m"
# Colormap (higher values = darker color). Examples: 'Reds', 'Oranges', 'Greys', 'Blues', 'YlOrRd'
COLORMAP = 'Reds'
# --- End Configuration ---

def plot_data(gdf, output_image_path, color_field=None, bounds=None):
    """
    Plots GeoDataFrame polygons on a map, optionally color-coded by a field,
    optionally filtered/zoomed by bounds, and saves the plot to a file.
    """

    if gdf.empty:
        print("Input GeoDataFrame is empty. Nothing to plot or save.")
        return

    print(f"Total features loaded: {len(gdf)}")
    gdf_to_plot = gdf
    map_extent = None # Used to set map view

    # --- Color Field Validation ---
    use_color_coding = False
    if color_field:
        if color_field not in gdf.columns:
            print(f"Warning: Specified color field '{color_field}' not found in GeoJSON properties. Plotting without color-coding.")
        else:
            # Attempt to convert column to numeric, coercing errors to NaN
            numeric_col = pd.to_numeric(gdf[color_field], errors='coerce')
            if numeric_col.isnull().all():
                 print(f"Warning: Color field '{color_field}' contains no valid numeric data. Plotting without color-coding.")
            elif not pd.api.types.is_numeric_dtype(numeric_col):
                 # This check might be redundant after to_numeric, but good for clarity
                 print(f"Warning: Color field '{color_field}' is not a numeric type. Plotting without color-coding.")
            else:
                 # If we reached here, the field exists and has at least some numeric data
                 use_color_coding = True
                 # Replace original column with numeric version (with NaNs for non-numeric)
                 # This helps geopandas handle the plotting correctly
                 gdf_to_plot[color_field] = numeric_col
                 print(f"Using field '{color_field}' for color-coding with colormap '{COLORMAP}'.")
                 # Report how many values are missing/non-numeric in the color field
                 missing_count = numeric_col.isnull().sum()
                 if missing_count > 0:
                     print(f"Note: {missing_count} features have missing or non-numeric values in '{color_field}' and will be plotted with the 'missing' style.")

    # --- End Color Field Validation ---


    if bounds:
        print(f"Filtering data view to bounds: {bounds}")
        min_lon, min_lat, max_lon, max_lat = bounds
        map_extent = bounds # Use bounds to set the map view
        # Optional filtering logic can go here if needed

    print(f"Generating plot for saving to: {output_image_path}")
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))

    # --- Plot the Polygons ---
    plot_kwargs = {
        'ax': ax,
        'edgecolor': 'grey', # Use a neutral edge color
        'linewidth': 0.5,
        'alpha': 0.75, # Adjust alpha as needed
    }

    if use_color_coding:
        plot_kwargs.update({
            'column': color_field,
            'cmap': COLORMAP,
            'legend': True,
            'legend_kwds': {'label': f"{color_field} Value",
                            'orientation': "horizontal",
                            'shrink': 0.6}, # Adjust legend size/position
            # Define how to style features with missing values in the color_field
            'missing_kwds': {
                "color": "lightgrey",
                "edgecolor": "darkgrey",
                "hatch": "//",
                "label": "Missing values",
            },
        })
        # Remove fixed facecolor if using color coding
        # plot_kwargs.pop('facecolor', None) # No need as it wasn't added yet
    else:
        # Default style if not color-coding
        plot_kwargs['facecolor'] = 'lightcoral'
        plot_kwargs['label'] = 'Geographic Areas' # Generic label

    gdf_to_plot.plot(**plot_kwargs)
    # --- End Polygon Plot ---


    # Set map extent if bounds were provided
    if map_extent:
        min_lon, min_lat, max_lon, max_lat = map_extent
        ax.set_xlim(min_lon, max_lon)
        ax.set_ylim(min_lat, max_lat)

    # Add basemap using contextily
    try:
        ctx.add_basemap(ax, crs=gdf_to_plot.crs.to_string(), source=ctx.providers.OpenStreetMap.Mapnik)
    except Exception as e:
         print(f"Warning: Could not add basemap. Error: {e}")


    # Customize plot
    plot_title = f'Crime Summary Areas by {color_field}' if use_color_coding else 'Crime Summary Areas'
    if bounds:
        plot_title += ' within Los Angeles Area Bounds'
    ax.set_title(plot_title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # Legend is handled by geopandas plot if color_coding is used
    # if not use_color_coding: plt.legend() # Add legend only if not color coding (if label was set)

    plt.tight_layout() # Adjust layout before saving

    # --- Save the plot instead of showing it ---
    try:
        plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
        print(f"Plot successfully saved to: {output_image_path}")
    except Exception as e:
        print(f"Error saving plot to {output_image_path}: {e}")

    # Close the plot figure to free memory
    plt.close(fig)
    # --- End Saving ---


if __name__ == "__main__":

    # python -m src.visualize_crime_data data/serious-violent-crimes_all.geojson
    # python -m src.visualize_crime_data data/serious-violent-crimes_all.geojson --no_la_filter

    parser = argparse.ArgumentParser(description="Visualize polygon data from GeoJSON, optionally color-coded by a numeric field, and save map as image.")
    parser.add_argument(
        "input_file",
        help="Path to the input GeoJSON file containing polygon data."
    )
    parser.add_argument(
        "--color_field",
        default=DEFAULT_COLOR_FIELD,
        help=f"Field name in GeoJSON properties to use for color-coding polygons. Defaults to '{DEFAULT_COLOR_FIELD}'."
    )
    parser.add_argument(
        "--no_la_filter",
        action="store_true",
        help="Plot all data in the file without zooming to LA bounds."
    )
    parser.add_argument(
        "--output_dir",
        default=IMG_OUTPUT_DIR,
        help=f"Directory to save the output image file. Defaults to '{IMG_OUTPUT_DIR}'"
    )


    args = parser.parse_args()

    # --- Prepare Output ---
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found at '{args.input_file}'")
        exit(1)

    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        except OSError as e:
            print(f"Error creating output directory {output_dir}: {e}")
            exit(1)

    base_name = os.path.basename(args.input_file)
    name_without_ext = os.path.splitext(base_name)[0]
    # Include color field name in output filename if not default
    color_suffix = f"_by_{args.color_field.lower()}" if args.color_field != DEFAULT_COLOR_FIELD else ""
    output_image_filename = f"{name_without_ext}{color_suffix}.png"
    output_image_path = os.path.join(output_dir, output_image_filename)
    # --- End Prepare Output ---


    # Load data
    try:
        print(f"Loading data from: {args.input_file}")
        data_gdf = gpd.read_file(args.input_file)
    except Exception as e:
        print(f"Error loading GeoJSON file: {e}")
        exit(1)

    # Determine bounds for plotting view
    plot_bounds = None if args.no_la_filter else LA_BOUNDS

    # Plot the data and save the image, passing the color field argument
    plot_data(data_gdf, output_image_path, color_field=args.color_field, bounds=plot_bounds)

    print("Visualization script finished.")