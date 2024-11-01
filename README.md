# Satellite Coverage Analysis

This repository provides tools to calculate unique coverage areas for a constellation of satellites. The project uses satellite Two-Line Element (TLE) data, geometric and geographic operations, and efficient spatial data processing methods to compute total and unique coverage areas. It leverages libraries like `pandas`, `geopandas`, `shapely`, and `pyorbital` for streamlined satellite and geographic data manipulation.

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

To get started with this project, ensure you have Python installed, along with the required packages outlined below. The repository is set up to handle TLE data, calculate satellite coverage, and assess spatial overlap.

## Features

1. **TLE Download and Conversion**: Download TLE data for satellite constellations, such as Starlink, directly from CelesTrak and convert it into a structured pandas DataFrame.
2. **Satellite Position Calculation**: Calculate real-time latitude, longitude, and altitude of satellites using TLE data.
3. **Coverage Area Calculation**: 
   - **Total Coverage Area**: Computes coverage radius and area based on latitude, longitude, and altitude.
   - **Unique Coverage Area**: Uses BallTree spatial indexing to determine overlapping satellites, calculating the union of all coverage areas to avoid double-counting.

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/satellite-coverage-analysis.git
    cd satellite-coverage-analysis
    ```

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