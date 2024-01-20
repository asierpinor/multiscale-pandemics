# USA data sources

Further information about data sources and manipulations can be found in Ref.~[PAPER]. See also `data_to_Tmatrix.md` and `usa_data_Tij.py`.

The specific data files used can be found in the `data` folder.

## Population data

We use 2021 population by county data from the US Census Bureau (https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html#par_textimage). Specifically, we used the file from https://www2.census.gov/programs-surveys/popest/datasets/2020-2021/counties/totals/. The data used can be found in `data/UScensus_data`.

## Commuting data

We use 2011-2015 residence-to-work commuting data from the US Census Bureau, see https://www.census.gov/data/tables/2015/demo/metro-micro/commuting-flows-2015.html. The data used can be found in `data/UScensus_data`.

## Air travel data

We use 2019 passenger boarding data by airport from the FAA and consider only enplanements in commercial airports, see https://www.faa.gov/airports/planning_capacity/passenger_allcargo_stats/passenger/previous_years#2019. The data used can be found in `data/FAA_data`.

## Cartographic boundary data

For figures we use 2022 cartographic boundary data (shapefile) from the US Census (https://www.census.gov/geographies/mapping-files/time-series/geo/cartographic-boundary.html), Canada states 2021 data from Statistics Canada (https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?Year=21), Mexico states data from OCHA (https://data.humdata.org/dataset/cod-ab-mex?), world countries data from Natural Earth (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-0-countries/), and continents data from ArcGIS Hub (https://hub.arcgis.com/datasets/esri::world-continents/about).
The data used can be found in `data/geodata`.

## City to county data

We use 2022 city-to-county data from Simple Maps, see US Cities Database in https://simplemaps.com/data/us-cities.
The data used can be found in `data/Simplemaps_data`.



