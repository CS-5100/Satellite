# Satellite Coverage Analysis

This repository provides tools to determine optimal satellite positions for achieving maximum area coverage on Earth. By using Two-Line Element (TLE) data, geographic transformations, and spatial data processing, this project analyzes the unique coverage area of satellite constellations. It provides functions to calculate each satellite’s position in 3D Earth-Centered coordinates and to compute coverage overlap, enabling analysis and optimization of satellite configurations for efficient global coverage.

## Table of Contents

- [Getting Started](#getting-started)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Download and Process TLE Data](#download-and-process-tle-data)
  - [Calculate Inter-Satellite Distance](#calculate-inter-satellite-distance)
  - [Compute Unique Coverage Area](#compute-unique-coverage-area)
- [Files](#files)
- [Contributing](#contributing)
- [License](#license)

---

## Getting Started

To get started with this project, ensure you have Python installed, along with the required packages outlined in the `requirements.txt` file. The repository is set up to handle TLE data, calculate satellite coverage, and assess spatial overlap.

## Features

1. **TLE Download and Conversion**: Download TLE data for satellite constellations, such as Starlink, directly from CelesTrak and convert it into a structured pandas DataFrame.
2. **Satellite Position Calculation**: Calculate real-time latitude, longitude, and altitude of satellites using TLE data.
3. **Coverage Area Calculation**: 
   - **Total Coverage Area**: Computes coverage radius and area based on latitude, longitude, and altitude.
   - **Unique Coverage Area**: Uses BallTree spatial indexing to determine overlapping satellites, calculating the union of all coverage areas to avoid double-counting.

## Installation

1. Clone the repository:
    ```bash
    cd /your/preferred/location/
    git clone https://github.com/CS-5100/Satellite.git
    cd Satellite/
    ```
    > **Note**: GitHub has removed [password-based authentication](https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls) when cloning with HTTPS URLs. Instead, use **GitHub CLI** or **personal access tokens (PAT)** for authentication when prompted.

2. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

   **Note:** The repository requires external libraries, including:
   - `pandas`
   - `numpy`
   - `requests`
   - `pyorbital`
   - `shapely`
   - `geopandas`
   - `scikit-learn`
   - `networkx`

## Usage

### Download and Process TLE Data

To download the latest TLE data and store it in a DataFrame:

```python
from tle_processing import download_tle_data_to_dataframe

url = "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=starlink&FORMAT=tle"
df = download_tle_data_to_dataframe(url)
print(df.head())
```

### Calculate Inter-Satellite Distance

Given two satellite positions, the calculate distance script calculates the Cartesian (Euclidean) distance between two satellites in Earth-Centered, Earth-Fixed (ECEF) coordinates.

```python
from satellite_distance import calculate_distance

# Example positions in Latitude, Longitude, and Altitude format
position1 = (550, 53, 0)    # Satellite 1's position
position2 = (540, 50, 20)   # Satellite 2's position

# Calculate the distance between two satellites
distance = calculate_distance(position1, position2)
```

### Computer Unique Coverage Area

To calculate the unique coverage area:

```python
from calculateUniqueCA_NEW import calculate_unique_coverage_area

# Assuming 'df' is your DataFrame with columns: 'Satellite Name', 'Latitude', 'Longitude', and 'Coverage Area (km^2)'
unique_coverage_area = calculate_unique_coverage_area(df)
print(f"Total unique coverage area: {unique_coverage_area} km²")
```

This function uses BallTree spatial indexing to compute unique coverage areas by grouping overlapping satellites and calculating the union of their coverage areas.