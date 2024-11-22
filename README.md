# Satellite Coverage Analysis

The repository contains a set of tools aimed at optimizing the placing of the Starlink constellation in order to achieve maximum area coverage. Using Two-Line Element data, geographic transformations, and spatial data processing, the project analyzes and optimizes the unique area covered by the satellite network.

## Table of Contents

- [Project Overview](#project-overview)
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
## Project Overview

The SpaceX Starlink network is currently the most extensive and technologically advanced satellite constellation in low earth orbit, offering high-speed, low-latency broadband internet worldwide. This work focuses on improving the performance of Starlink by solving the optimization problem of satellite positioning to achieve maximum area coverage.

## Getting Started

To get started with this project, ensure you have Python installed.

### Python Installation
This project requires Python 3.9 or 3.10. You can download and install Python for:

Windows: Visit the [Python](https://www.python.org/downloads/) website and download the installer.

macOS: Using Homebrew: 
```bash
brew install python@3.9
```
Linux: Most distributions come with Python pre-installed. If not, use your package manager: 
```bash
# Ubuntu/Debian:
sudo apt-get install python3.9
```
```bash
# Fedora: 
sudo dnf install python3.9
```

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
    The required librabries to preprocess the TLE data and run the local search algorithms would be installed. 
    
    Note: It is recommended to use a virtual environment to install the libraries to prevent any discrepancies with the existing libraries on your system.

## Usage

### Download and Process TLE Data

To download the latest TLE data and store it as a GeoDataFrame:

```python
import tle_processing as tlp

# Download and convert TLE data to GeoDataFrame for existing satellites
current_tle = tlp.download_current_tles_as_list()
tle_df = tlp.tles_to_dataframe(raw_tle_list=current_tle)
satellite_gdf = tlp.tle_dataframe_to_geodataframe(tle_df)
```
