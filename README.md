# Snow-Crab-Condition
Data and analyses for NPRB funded project #1911 "Bering Sea snow crab lipid condition metrics in relation to temperature and recruitment dynamics". This dataset includes immature snow crab sampled for hepatopancreas on 2019, 2021, 2022 and 2023 EBS and NBS NOAA bottom trawl surveys. Fatty acid analyses from hepatopancreas samples were used to explore spatiotemporal variation in condition and lipid biomarkers in relation to temperature. 

# Objectives:
1)	Determine spatial and interannual variation in body condition (i.e. total lipid content) of juvenile snow crab in the Bering Sea 
2)	Relate annual variation in lipid condition metrics to Bering Sea temperatures and annual estimates of snow crab recruitment in order to understand the importance of an energetic-based mechanism in linking thermal conditions to recruitment
3)	Assess the performance of rapid condition metric in estimating energetic condition of field-collected snow crab 

# METADATA:
================================================

Bibliographic
-------------

| Published     | 02/02/2023   |
| ------------- | ------------- |
| Keywords      | snow crab |
|               |   Bering Sea |
|               |  fatty acids |
|               | condition |

 

Coverage
--------

### Temporal

| Begin    | 2019-06-01 |
| ------------- | ------|
| End   | 2023-08-31 |

 
 

### Spatial

| LME     |                     |
| ------------- | ------------- |
|                | Eastern Bering Sea |
|               |   Northern Bering Sea |

![Rplot]([https://user-images.githubusercontent.com/59858752/216689865-6766501b-98e4-420e-98cd-40a75e7a4ca1.png](https://github.com/Erin-Fedewa-NOAA/Snow-Crab-Condition/blob/main/figures/data%20exploration/n_year.png?raw=true))


Attributes
----------
Two master datasets have been produced for further analyses via the "append haul" script. 1) total_FA_master.csv includes all snow crab biometric, total fatty acid and haul level data, with data attributes listed below, and 2) FA_biomarker_master.csv includes individual fatty acid biomarkers for each hepatopancreas sample by % weight, per wet weight (micrograms/mg), and per dry weight (micrograms/mg). Haul and biometric data in the attribute table below are included as well.

| Name    |    Description   |   Unit    |
| ------- | ---------------- | ---------- |
| `cruise` | Cruise ID for Bering Sea bottom trawl surveys. YYYY-01 indicates EBS surveys, YYYY-02 indicates NBS surveys. See RACEBASE or AKFIN for additional NOAA cruise metadata  |   ID Code
| `gis_station`   | Alpha-numeric designation for the station established in the design of AFSC standardized surveys | numeric/text                                      
| `area_swept`   |   Unit of effort for AFSC bottom trawl surveys: computed by distance towed*mean net width   | numeric, in ha
|  `cpue`   |   Station-level snow crab density, calculated as CPUE"   |   numeric, crab/nmi^2
|  `vial_id`    |    Unique ID for hepatopancreas sample by year. AKK=Alaska Knight samples, V=Vesteraalen samples   |   numeric/text
|  `vessel`  |     ID number of the vessel used to collect data for that haul associated with vessel name    |   numeric
|  `haul`      |  Uniquely identifies a sampling event (haul) within an AFSC cruise. It is a sequential number, in chronological order of occurrence |  numeric
|  `sex` | Sex of specimen sampled. 1=Male, 2-Female   |  numeric
|  `cw`  |  Carapace width of specimen sampled | numeric, in mm
|  `ch_cc`   |   Chela height (males) or clutch code (females) used to determine maturity. 0 clutch code for females signifies 000, or immature   | numeric, in mm
 | `crab_wgt`  |   Whole crab weight of specimen sampled (prior to dissection). Only individuals with no missing limbs weighed | numeric, in g
|  `hepato_wwt`    | Total wet weight of hepatopancreas (subtracted from weight of weigh boat)  | numeric, in g
|  `hepato_dwt`    |Total dry weight of hepatopancreas after drying at 70C for 48hrs (subtracted from weight of weigh boat)  | numeric, in g
|  `maturity`  |  Maturity of specimen sampled. 0=Immature, 1=Mature"  | numeric                                                                                        
| 'mid_latitude'       |   Latitude of specimen collection. Designates latitude at start of haul for AFSC standardized surveys    | numeric
|  `mid_longitude`    | Longitude of specimen collection. Designates longitude at start of haul for AFSC standardized surveys | decimal degree
 | `bottom_depth`    |    Bottom depth at station for AFSC standardized surveys  | numeric, in m
|  `gear_temperature`   |    Bottom temperature at sampling station | degree C
 | `year`     |        Year of specimen collection | numeric
|  `sample_region`        |   Pre-defined spatial sampling strata in the Bering Sea for hepatopancreas collections. Regions 1-6 correspond with EBS BSIERP regions, Regions 7-9 correspond with NBS BSIERP regions  | numeric
|  `lme` |   Large marine ecosystem: EBS=eastern Bering Sea, NBS=northern Bering Sea   |  text
|  `DWT_WWT` |   hepatopancreas dry weight/hepatopancreas wet weight   |  numeric, in g
|  `Perc_DWT` |   hepatopancreas dry weight/ hepatopancreas wet weight x 100; percentage hepatopancreas dry weight    |  percentage
|  `Total_FA_Conc_WWT` |   total fatty acid wet weight concentration of hepatopancreas sample   |  mg FA/ g WWT
|  `Total_FA_Conc_DWT` |   total fatty acid dry weight concentration of hepatopancreas sample   |  mg FA/ g DWT
|  `Total_FA` |   Total fatty acid concentration/hepatopancreas dry weight:wet weight ratio   |  mg/g DWT
|  `WWT_DWT` |   hepatopancreas wet weight/dry weight    |  numeric, in g

:::

Distribution
------------

  File                                              Format    
  ------------------------------------------------- -------- -------------------------------------------------
  `total_FA_master`   `csv`    [Download](https://github.com/Erin-Fedewa-NOAA)
