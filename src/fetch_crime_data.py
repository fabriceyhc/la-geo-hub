# fetch_crime_data.py
import requests
import json
import os
import argparse
from urllib.parse import urlparse
from datetime import datetime
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, Polygon

# --- Configuration ---
BASE_URL = "https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/Serious_Violent_Crimes/FeatureServer/0/query"
# !!! IMPORTANT: Verify this is the correct date field name IF using date filters !!!
# This field MUST exist in the specific layer you query if you provide date arguments.
# Based on your previous message, Layer 0 does NOT have a date field.
DATE_FIELD = "DateReported" # <--- Placeholder - CHANGE IF NEEDED FOR A LAYER WITH DATES
OUTPUT_DIR = "data"
# ArcGIS often uses timestamps. We'll append time for broader compatibility.
DATE_FORMAT_ARCGIS = "%Y-%m-%d %H:%M:%S"
# --- End Configuration ---

def fetch_data(start_dt=None, end_dt=None):
    """
    Fetches data from the ArcGIS service.
    If start_dt and end_dt are provided, filters by date range.
    Otherwise, fetches all data (where=1=1).
    """
    params = {
        'outFields': '*',
        'outSR': '4326',  # WGS 84 Latitude/Longitude
        'f': 'json',
        # Add resultOffset and resultRecordCount here if pagination is needed
        # 'resultRecordCount': 2000 # Example limit per request (max often 2000 for ArcGIS Online)
    }

    # --- Construct WHERE clause ---
    if start_dt and end_dt:
        # Only construct date filter if BOTH dates are provided
        print(f"--- Date filtering enabled for field: {DATE_FIELD} ---")
        # Format dates for ArcGIS WHERE clause (using timestamps)
        start_str = start_dt.strftime(DATE_FORMAT_ARCGIS)
        # Add time to end_dt to include the whole day
        end_dt_inclusive = datetime.combine(end_dt.date(), datetime.max.time())
        end_str = end_dt_inclusive.strftime(DATE_FORMAT_ARCGIS)

        # Construct the WHERE clause for date filtering using 'timestamp' keyword
        # Example: "IncidentDate >= timestamp '2025-03-29 00:00:00' AND IncidentDate <= timestamp '2025-04-28 23:59:59'"
        where_clause = f"{DATE_FIELD} >= timestamp '{start_str}' AND {DATE_FIELD} <= timestamp '{end_str}'"
        print(f"Using WHERE clause: {where_clause}")
        params['where'] = where_clause
    else:
        # Fetch all data if dates are not provided
        print("--- No date filters provided. Fetching all data. ---")
        where_clause = "1=1"
        params['where'] = where_clause
    # --- End WHERE clause construction ---


    all_features = []
    offset = 0
    exceeded_limit = True
    record_limit = params.get('resultRecordCount', 2000) # Use specified limit or default guess

    print(f"Fetching data from {BASE_URL}...")

    # Pagination handling loop
    while exceeded_limit:
        if offset > 0:
            params['resultOffset'] = offset
            print(f"Fetching next batch (offset={offset})...")

        try:
            response = requests.get(BASE_URL, params=params, timeout=90)
            response.raise_for_status()
            data = response.json()

        except requests.exceptions.RequestException as e:
            print(f"Error fetching data: {e}")
            return None
        except json.JSONDecodeError:
            print("Error decoding JSON response.")
            print("Response text:", response.text[:500])
            return None

        # Check for errors reported within the JSON response
        if data.get('error'):
            # Specific check for invalid date field when filtering is attempted
            if start_dt and end_dt and 'invalid field' in str(data['error']).lower() and DATE_FIELD in str(data['error']):
                 print(f"\nAPI Error: The specified DATE_FIELD ('{DATE_FIELD}') is invalid for this layer.")
                 print("Cannot perform date filtering. Check the layer's fields.")
                 print("If you want all data, run the script without --start_date and --end_date arguments.\n")
            else:
                print(f"API Error: {data['error'].get('message', 'Unknown error')}")
                print(f"Details: {data['error'].get('details', [])}")
            return None # Stop fetching on error

        features = data.get('features', [])
        if features:
            all_features.extend(features)
            print(f"Fetched {len(features)} features in this batch.")
        else:
            # If offset is 0 and no features, it means no data for the query
            if offset == 0:
                 print("No features found for the specified query.")
            else:
                 print("No more features found in this batch.")
            break # Exit loop if no features returned

        # Check if the transfer limit was exceeded OR if we received exactly the limit (common case)
        exceeded_limit = data.get('exceededTransferLimit', False) or (len(features) == record_limit)

        if exceeded_limit:
            offset += len(features)
            if not features: # Safety break if exceededLimit is true but no features returned
                 print("Warning: Pagination indicated more data, but no features were returned. Stopping.")
                 break
        else:
             print(f"Finished fetching. Total features: {len(all_features)}")
             break # Exit loop if limit wasn't exceeded

    return all_features


def save_features_to_geojson(features, filename):
    """
    Converts features (Points OR Polygons) to GeoDataFrame and saves as GeoJSON.
    """
    if not features:
        print("No features to save.")
        return False

    print(f"Processing {len(features)} features for saving...")
    try:
        # Extract attributes - Handle potential missing attributes key
        attributes = [f.get('attributes', {}) for f in features]

        # --- Geometry Extraction (Handles Points and Polygons) ---
        geometry = []
        valid_feature_indices = []

        # Determine geometry type from the first feature if possible (assume consistent type)
        # This isn't strictly necessary but can provide context
        geom_type = None
        if features and features[0].get('geometry'):
            if 'rings' in features[0]['geometry']:
                geom_type = 'Polygon'
            elif 'x' in features[0]['geometry'] and 'y' in features[0]['geometry']:
                geom_type = 'Point'
        print(f"Attempting to process geometry type: {geom_type or 'Unknown/Mixed'}")


        for i, feature in enumerate(features):
            geom_data = feature.get('geometry')
            processed = False
            if geom_data:
                # Try processing as Polygon first
                if 'rings' in geom_data and geom_data['rings']:
                    try:
                        # Esri rings format: [[ [x1,y1], [x2,y2], ... ], [ [hole1_x1,y1], ... ]]
                        # Shapely Polygon takes shell (outer ring) then optional holes (inner rings)
                        shell = geom_data['rings'][0]
                        holes = geom_data['rings'][1:] if len(geom_data['rings']) > 1 else None
                        # Basic validation: Need at least 3 unique points for a valid ring
                        if len(shell) >= 4 and shell[0] == shell[-1]: # Esri rings are usually closed
                             poly = Polygon(shell, holes=holes)
                             if poly.is_valid:
                                 geometry.append(poly)
                                 valid_feature_indices.append(i)
                                 processed = True
                             else:
                                 print(f"Warning: Skipping feature {i} due to invalid polygon geometry after creation.")
                        # else: print(f"Debug: Skipping feature {i} polygon - insufficient points or not closed.")
                    except Exception as e:
                        print(f"Warning: Skipping feature {i} polygon due to error during creation: {e}")

                # If not processed as Polygon, try processing as Point
                elif not processed and 'x' in geom_data and 'y' in geom_data:
                     if geom_data['x'] is not None and geom_data['y'] is not None:
                         try:
                             pt = Point(geom_data['x'], geom_data['y'])
                             geometry.append(pt)
                             valid_feature_indices.append(i)
                             processed = True
                         except Exception as e:
                             print(f"Warning: Skipping feature {i} point due to error during creation: {e}")
                     # else: print(f"Debug: Skipping feature {i} point - null coordinates.")

            # if not processed: print(f"Debug: Skipping feature {i} - No recognizable geometry found.")


        if not geometry:
             print("Error: No valid geometries (Point or Polygon) could be extracted from the fetched features.")
             return False

        # Ensure attributes DataFrame matches the valid geometries
        if len(valid_feature_indices) != len(features):
             print(f"Warning: {len(features) - len(valid_feature_indices)} features lacked valid/supported geometry and were excluded.")
             # Filter attributes to only include those with valid geometry
             attributes_df = pd.DataFrame([attributes[i] for i in valid_feature_indices])
        else:
             attributes_df = pd.DataFrame(attributes)

        # Check if attributes_df is empty after filtering
        if attributes_df.empty and geometry:
             print("Error: No attributes found for the valid geometries. Cannot create GeoDataFrame.")
             return False
        elif not geometry and not attributes_df.empty:
             print("Error: No geometries found for the attributes. Cannot create GeoDataFrame.")
             return False
        elif attributes_df.empty and not geometry:
             print("Error: No attributes or geometries processed.") # Should be caught earlier
             return False

        # --- End Geometry Extraction ---

        gdf = gpd.GeoDataFrame(attributes_df, geometry=geometry, crs="EPSG:4326")
        print(f"Created GeoDataFrame with {len(gdf)} features.")

        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        # Save to GeoJSON
        gdf.to_file(filename, driver='GeoJSON')
        print(f"Data successfully saved to: {filename}")
        return True

    except Exception as e:
        # Catch broader errors during processing/saving
        import traceback
        print(f"An unexpected error occurred in save_features_to_geojson: {e}")
        print(traceback.format_exc()) # Print detailed traceback
        return False


if __name__ == "__main__":

    # python -m src.fetch_crime_data 
    # python -m src.fetch_crime_data --layer_url https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/Serious_Violent_Crimes/FeatureServer/1/query
    # python -m src.fetch_crime_data --start_date 2025-04-01 --end_date 2025-04-15 --layer_url https://services.arcgis.com/RmCCgQtiZLDCtblq/arcgis/rest/services/Serious_Violent_Crimes/FeatureServer/1/query

    
    parser = argparse.ArgumentParser(description="Fetch crime data from ArcGIS REST API. Optionally filter by date range (if supported by layer and DATE_FIELD is correct). Saves as GeoJSON.")

    parser.add_argument(
        "--start_date",
        default=None, # Default is None unless user provides it
        help="Start date in YYYY-MM-DD format (inclusive). Requires --end_date."
    )
    parser.add_argument(
        "--end_date",
        default=None, # Default is None unless user provides it
        help="End date in YYYY-MM-DD format (inclusive). Requires --start_date."
    )
    parser.add_argument(
        "--output_dir",
        default=OUTPUT_DIR,
        help=f"Directory to save the output GeoJSON file. Defaults to '{OUTPUT_DIR}'"
    )
    parser.add_argument(
        "--layer_url",
        default=BASE_URL,
        help="Full URL to the Feature Service Layer Query endpoint (optional)."
    )


    args = parser.parse_args()

    # Update BASE_URL if provided via argument
    if args.layer_url != BASE_URL:
        BASE_URL = args.layer_url
        print(f"Using custom layer URL: {BASE_URL}")


    start_datetime = None
    end_datetime = None

    # --- Validate Date Arguments ---
    if args.start_date and not args.end_date:
        parser.error("--start_date requires --end_date.")
    if args.end_date and not args.start_date:
        parser.error("--end_date requires --start_date.")

    if args.start_date and args.end_date:
        try:
            # Convert string dates to datetime objects (start of day)
            start_datetime = datetime.strptime(args.start_date, '%Y-%m-%d')
            end_datetime = datetime.strptime(args.end_date, '%Y-%m-%d')
        except ValueError:
            parser.error("Invalid date format. Please use YYYY-MM-DD.")

        if start_datetime > end_datetime:
            parser.error("Start date cannot be after end date.")
    # --- End Date Validation ---


    # Fetch data - pass potentially None dates
    fetched_features = fetch_data(start_datetime, end_datetime)

    # Save data if fetched successfully
    if fetched_features:
# --- Determine filename ---
        try:
            # Parse the URL used for fetching (BASE_URL reflects args.layer_url if provided)
            parsed_url = urlparse(BASE_URL)
            # Split path into components, removing empty strings (like leading '/')
            path_parts = [part for part in parsed_url.path.split('/') if part]

            # Find the index of 'services'
            services_index = path_parts.index('services')

            # The service name should be the part immediately after 'services'
            if services_index + 1 < len(path_parts):
                service_name_raw = path_parts[services_index + 1]
                # Basic cleaning: lower case, replace spaces/underscores with hyphens
                service_name = service_name_raw.lower().replace('_', '-').replace(' ', '-')
            else:
                # Fallback if path ends right after 'services'
                service_name = 'features'
                print(f"Warning: Could not determine service name component after 'services' in URL path. Using fallback '{service_name}'.")

        except (ValueError, IndexError, AttributeError):
            # General fallback if URL parsing or 'services' finding fails
            service_name = 'features'
            print(f"Warning: Could not automatically determine service name from URL '{BASE_URL}'. Using fallback '{service_name}'.")


        # Now create the filename using the determined service_name
        if start_datetime and end_datetime:
            # Create filename with dates
            start_str_file = start_datetime.strftime('%Y%m%d')
            end_str_file = end_datetime.strftime('%Y%m%d')
            output_filename = os.path.join(args.output_dir, f"{service_name}_{start_str_file}_{end_str_file}.geojson")
        else:
            # Create generic filename for all data
            output_filename = os.path.join(args.output_dir, f"{service_name}_all.geojson")
        # --- End Filename Determination ---

        # (The call to save_features_to_geojson remains the same)
        save_features_to_geojson(fetched_features, output_filename)
    else:
        print("\nScript finished: Could not fetch or process data.")