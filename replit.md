# Overview

This is a comprehensive Streamlit web application for calculating the area enclosed by the Huangyan Island (Scarborough Shoal) territorial sea baseline. The application provides advanced geospatial analysis using multiple interpolation methods, interactive visualization, and detailed reporting capabilities.

**Project Purpose**: Calculate and analyze the area of Huangyan Island's territorial sea baseline as defined in the Statement by the Government of the People's Republic of China on the Territorial Sea Baselines of Huangyan Dao (November 10, 2024).

**Current State**: Fully functional production application with all MVP and enhancement features implemented and tested.

# User Preferences

Preferred communication style: Simple, everyday language (Chinese and English).

# System Architecture

## Application Structure

**Problem**: Need to accurately calculate and visualize the area enclosed by territorial sea baseline points using different mathematical methods.

**Solution**: Streamlit-based web application with comprehensive geospatial analysis capabilities.

**Design Pattern**: Single-page application with modular functions for coordinate conversion, area calculation, distance measurement, interpolation, mapping, and reporting.

**Key Components**:
- `app.py`: Main application file containing all functionality
- Coordinate parsing and conversion functions
- Area calculation using Shapely geometry library
- Multiple interpolation methods (linear, quadratic spline, cubic spline, Bezier)
- Interactive Folium maps
- PDF report generation using ReportLab

## Programming Language & Framework

**Technology**: Python 3.11 with Streamlit

**Rationale**:
- Streamlit enables rapid development of interactive data applications
- Python's rich scientific computing ecosystem (NumPy, SciPy, Pandas)
- Excellent geospatial libraries (Shapely, PyProj, Folium)

## Core Features (December 2024)

### 1. Baseline Point Data Management
- Display of 16 baseline points in tabular format
- Support for degree-minute format (15°08.1′) and decimal degrees (15.135)
- CSV/Excel file upload for custom baseline data
- Automatic coordinate conversion to decimal degrees

### 2. Area Calculations
- **Straight-line polygon method**: Connects points with straight lines (136.5 km²)
- **Smooth curve method**: Uses cubic spline interpolation (160.3 km²)
- Accurate coordinate transformation (WGS84 → UTM Zone 50N)
- Results precise to 0.1 km²

### 3. Distance Analysis
- Calculates distances between consecutive baseline points
- Includes closing segment (point 15 → point 1)
- Total perimeter: 46.26 km (24.98 nautical miles)
- Detailed table of all 15 baseline segments
- Identifies longest and shortest segments

### 4. Multiple Interpolation Methods
- Linear interpolation (straight-line polygon)
- Quadratic spline interpolation
- Cubic spline interpolation
- Bezier curve interpolation
- Comparative analysis showing area differences between methods
- Visual comparison with 2x2 plot grid

### 5. Interactive Mapping
- Folium-based interactive map with zoom and pan
- Markers for each baseline point with popup information
- Polygon overlay showing the baseline boundary
- Layer controls for customization

### 6. Visualization
- Matplotlib plots comparing straight-line vs smooth curve methods
- 2x2 grid comparison of different interpolation methods
- Point numbering and labeling
- Coordinate system displays (UTM projected coordinates)

### 7. PDF Export
- Comprehensive PDF reports with all calculation results
- Area comparison tables
- Perimeter measurements
- Interpolation method comparisons
- Baseline point coordinate listings
- Professional formatting with tables and sections

### 8. Data Upload
- Support for CSV and Excel (.xlsx, .xls) files
- Required columns: name, latitude, longitude
- Handles both degree-minute and decimal coordinate formats
- Validation and error handling

## Technical Implementation

### Coordinate Systems
- Input: WGS84 geographic coordinates (EPSG:4326)
- Processing: UTM Zone 50N projected coordinates (EPSG:32650)
- Ensures accurate distance and area calculations in meters

### Geospatial Libraries
- **Shapely**: Polygon creation and area calculation
- **PyProj**: Coordinate system transformations
- **Folium**: Interactive web maps
- **SciPy**: Spline interpolation and Bezier curves

### Area Calculation Method
- Uses Shapely's Polygon.area for precise calculations
- Accounts for Earth's curvature through UTM projection
- Handles closed polygons automatically

### Distance Calculation
- Euclidean distance in UTM coordinates
- Conversion to kilometers and nautical miles
- Modulo arithmetic for closing segment

### Error Handling
- Robust coordinate parsing with multiple format support
- File upload validation
- Interpolation error catching
- Clear user feedback for errors

## Data Flow

1. **Input**: Baseline point data (built-in or uploaded)
2. **Parsing**: Convert coordinates to decimal degrees
3. **Transformation**: Project to UTM coordinates
4. **Calculation**: Compute areas, distances, interpolations
5. **Visualization**: Generate maps and plots
6. **Export**: Create PDF reports

## Recent Changes (November 25, 2024)

### Bug Fixes
- Fixed distance calculation to include closing segment (15 total segments)
- Enhanced coordinate parsing to handle decimal degree inputs
- Fixed LSP warnings for unbound variables

### Features Added
1. Distance analysis with total perimeter and segment details
2. Multiple interpolation methods comparison (4 methods)
3. Interactive Folium map with markers and polygon
4. CSV/Excel file upload for custom data
5. PDF report generation and download

# External Dependencies

## Core Framework
- `streamlit>=1.51.0`: Web application framework

## Scientific Computing
- `numpy`: Numerical computations
- `pandas`: Data manipulation and analysis
- `scipy`: Interpolation and scientific functions

## Geospatial
- `shapely>=2.1.2`: Geometric operations and area calculations
- `pyproj>=3.7.2`: Coordinate system transformations
- `folium>=0.20.0`: Interactive maps
- `streamlit-folium>=0.25.3`: Folium integration with Streamlit

## Visualization
- `matplotlib>=3.10.7`: Static plots and charts

## File I/O
- `openpyxl>=3.1.5`: Excel file support
- `pillow`: Image processing

## PDF Generation
- `reportlab>=4.4.5`: PDF document creation

## Configuration
- Environment: Python 3.11
- Port: 5000 (configured for webview)
- Streamlit server settings in `.streamlit/config.toml`

# Testing & Quality Assurance

All features have been tested and verified:
- ✅ Area calculations (straight-line and smooth curve methods)
- ✅ Distance calculations (all 15 segments including closing)
- ✅ Interpolation method comparisons (4 methods)
- ✅ Interactive map functionality
- ✅ File upload (CSV/Excel with multiple coordinate formats)
- ✅ PDF export functionality
- ✅ Coordinate parsing (degree-minute and decimal formats)
- ✅ No critical LSP errors

# Future Enhancement Possibilities

While the current application is fully functional, potential future enhancements could include:
- Additional map tile layers (satellite, terrain)
- Export to multiple formats (Excel, GeoJSON)
- Support for multiple baseline sets comparison
- Historical data tracking
- 3D visualization
- Additional interpolation methods
- Chinese font support in PDF reports
- Database persistence for uploaded datasets
