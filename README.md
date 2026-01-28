# Oran-Marine-Pollution-GEE

## Project Overview
This project presents a coastal environmental monitoring framework for the Oran coastline (Algeria) using Google Earth Engine.
It integrates optical, radar, and environmental data to detect and analyze marine pollution and coastal dynamics.

## Objectives
- Detect potential oil spill events using SAR data
- Identify algae blooms using spectral indices
- Analyze water turbidity in coastal zones
- Monitor environmental conditions near strategic ports

## Study Area
The study focuses on the coastal and marine area of the Wilaya of Oran, including a 20 km offshore buffer.

## Data Sources
- Sentinel-1 (SAR): oil spill detection
- Sentinel-2 (optical): algae and turbidity analysis
- ERA5: wind speed and direction
- ETOPO1: bathymetry

## Methodology
The workflow combines:
- Marine water masking
- Spectral index computation (NDWI, NDCI, FAI)
- SAR backscatter analysis with wind-based filtering
- Multi-source data integration for coastal monitoring

## Outputs
- Oil spill detection maps
- Algae bloom maps
- Turbidity maps
- Sea surface temperature visualization

## Tools & Technologies
- Google Earth Engine
- JavaScript
- Remote Sensing & GIS

## Repository Structure
